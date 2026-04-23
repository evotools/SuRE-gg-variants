#!/usr/bin/env Rscript

start_time <- Sys.time()

suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(dtplyr))

# Stop on any error immediately
options(error = function() {
  cat("Error occurred at:", Sys.time(), "\n")
  cat("Error message:", geterrmessage(), "\n")
  traceback()
  quit(status = 1)  # Exit with error code
})
#options(warn = 2)  # Turn warnings into errors

# Variants with >2 alleles can't be tested with the model used for the main data
# set, and need to be tested differently
# This script processes the fragment data from non-biallelic variants 
# (libraries combined in combine_libs.R and these variants split off)
# and makes the data chunks for the models
# Similar steps to prep_data.nf.R but some differences needed
# All chromosomes are processed together, as there aren't enough variants to make batch processing worthwhile

# Some variants are non-biallelic in the vcf (i.e. >1 alt allele) but actually only have
# two alleles present in the SuRE data (e.g. if ref allele isn't present in either sample)
# These are output as a separate chunk to go back into the main model testing 

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
out_pref <- args[2]
source(config_file)

# Set threads for data.table as a parameter given to the script
n_threads <- as.integer(args[3])
setDTthreads(n_threads)

message(sprintf("Before loading data: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Load in data from all chr, adding chromosome location to rows and genotype (assumes all GT files are staged)
multi_files <- list.files(pattern = 'multi.comb.rds')
load_rdata <- function(file) {
  onechr.multi <- readRDS(file)
  chrid <- str_remove(file, fixed(paste0(".", out_pref, ".multi.comb.rds")))
  onechr.multi <- onechr.multi %>% mutate(chr = chrid)
  gts_loc <- paste(chrid, out_pref, "var_details.tsv", sep = ".")
  gts <- read_tsv(gts_loc, col_types = "cciccc") %>% select(-Chr, -RefLength, -AltLengths) %>% pivot_longer(cols=-VarID, names_to = "library", values_to = "genotype")
  onechr.multi <- onechr.multi %>% left_join(gts, join_by(VarID, library))
}
multi_data <- list_rbind(map(multi_files, load_rdata))

message(sprintf("Loaded data, starting to calculate filtering stats: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Load in thresholds calculated across all chr for the biallelic variants
# # Keep genotype 0/0, use the values for 0/1 for all het sites and 1/1 for all homozygous alt (of any alt allele)
ratioThresh <- read_tsv(paste0(out_pref,".thresholds.ratio.tsv"), show_col_types = F)
covThresh <- read_tsv(paste0(out_pref,".thresholds.iPCR.tsv"), show_col_types = F)
impThresh <- read_tsv(paste0(out_pref,".thresholds.imp.tsv"), show_col_types = F)
thresh <- left_join(ratioThresh, covThresh) %>% left_join(impThresh) %>% 
  mutate(gt.label=if_else(genotype=="0/0", '0/0', if_else(genotype=="1/1", 'homalt', 'het')), .keep='unused', .after='library')

# Calculate various metrics for filtering:
# Total coverage and proportion imputed as for biallelic variants (0 if no imputed)
# Instead of ref proportion get the proportion of the first allele in the genotype (i.e. with lowest ID, since gts have been unphased).
# For homozygous alt then take 1- this to compare with the thresholds for 1/1.
# Find the total (non-imputed) input counts for the allele from hap1, the allele from hap2 if different, and the total for any alleles that don't appear in the genotype.
# Then calculate fraction of input with allele 1. (For hom variants only hap1 will get a match, hap2 will always be NA)
# Also caclulate the 'ok proportion' - the fraction of input with alleles that appear in the genotype (1-don't appear/total)

allcounts <- multi_data %>% ungroup() %>% group_by(VarID, VarAllele, library, genotype) %>% 
  summarise(total=sum(input))
from_seq <- multi_data %>% filter(VarCav == 0 | VarCav == "E")
fs_counts <- from_seq %>% ungroup %>% group_by(VarID, VarAllele, library, genotype) %>% 
  summarise(nonI=sum(input))
setDT(allcounts)
setDT(fs_counts)
counts <- fs_counts[allcounts, on = c("VarID", "library", "VarAllele", "genotype"), nomatch = 0]
counts[is.na(nonI), nonI := 0]
counts[, c("hap1", "hap2") := tstrsplit(genotype, "/", fixed = TRUE)]
counts[, which.allele := fcase(
  VarAllele == hap1, "all1",
  VarAllele == hap2, "all2",
  default = "other"
)]
counts[, gt.label := fcase(
  hap1 != hap2, "het",
  hap1 == "0", "0/0",
  default = "homalt"
)]

by.var_library <- dcast(
  counts,
  VarID + library + genotype + gt.label ~ which.allele,
  value.var = c("total", "nonI"),
  fun.aggregate = sum,
  fill = 0
)
by.var_library[, total := total_all1 + total_other + total_all2]
by.var_library[, totalNonI := nonI_all1 + nonI_all2 + nonI_other]
by.var_library[, refPropEquiv := ifelse(gt.label == "homalt", 1 - nonI_all1/totalNonI, nonI_all1/totalNonI)]
by.var_library[, okprop := 1 - nonI_other/totalNonI]
by.var_library[, impPropvar := 1 - totalNonI/total]

# Calculate imputed proportion by allele instead of by variant, and for hets add whichever allele has the higher proportion, for filtering
# Ignore non-genotype fragments (which means a few variants might not have any data if only non-gts)
counts.imp.byal <- counts %>% filter(which.allele == 'all1' | which.allele == 'all2') %>% 
  mutate(impProp = 1-nonI/total)
counts.imp.var <- counts.imp.byal %>% 
  pivot_wider(id_cols = c(VarID, library, genotype, gt.label), names_from = which.allele, values_from = impProp, names_prefix='impProp_') %>% 
  replace_na(list(impProp_all1=0, impProp_all2=0)) %>% 
  rowwise() %>% 
  mutate(maxImp = max(impProp_all1, impProp_all2), .keep='unused')
by.var_library <- by.var_library %>% left_join(counts.imp.var)

# add ok proportion to thresholds
# Also as a sense check calculate what the other cutoffs would be if done from these data 
# (i.e. check if very different from the biallelic data set)
multi.thresh <- by.var_library %>% ungroup %>% group_by(library, gt.label) %>% 
  summarise(meanRatio=mean(refPropEquiv, na.rm=TRUE), downRatio=meanRatio-2.58*sd(refPropEquiv, na.rm=TRUE), upRatio=meanRatio+2.58*sd(refPropEquiv, na.rm=TRUE), meanImp=mean(impPropvar, na.rm=TRUE), upImp=meanImp+2.58*sd(impPropvar, na.rm=TRUE), meanOk=mean(okprop, na.rm=TRUE), downOk=meanOk-2.58*sd(okprop, na.rm=TRUE))
multi.cov <- by.var_library %>% ungroup %>% group_by(library, gt.label) %>% 
  filter(total < 1000) %>% summarise (meanCov=mean(total, na.rm=TRUE), downCov=meanCov-2.58*sd(total, na.rm=TRUE), upCov = meanCov+2.58*sd(total, na.rm=TRUE))
thresh <- thresh %>% left_join(multi.thresh %>% select(library, gt.label, meanOk, downOk))

message(sprintf("Calculated filtering stats, starting to do plots: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Plot distributions with the thresholds used (from the biallelic data) and those calculated from these data
# For imputed plot a value for each allele (but on same plot)
plot1 <- ggplot(counts.imp.byal, mapping=aes(x=impProp)) + geom_density() + 
  facet_wrap(vars(library, gt.label), scales='free') +
  geom_vline(data = thresh, aes(xintercept = meanImp), linetype='dashed') +
  geom_vline(data = thresh, aes(xintercept = upImp), linetype='dashed', colour='red') +
  geom_vline(data = multi.thresh, aes(xintercept = meanImp), linetype='dotted') + 
  geom_vline(data = multi.thresh, aes(xintercept = upImp), linetype='dotted', colour='pink') +
  labs(title = "Proportion of fragments that are imputed", subtitle="Dashed=from main data (used for filtering); dotted=from multi vars for comparison\nBlack=mean, red/blue=upper/lower thresholds")

plot2 <- ggplot(by.var_library, mapping=aes(x=refPropEquiv)) + geom_density() + 
  facet_wrap(vars(library, gt.label), scales='free') +
  geom_vline(data = thresh, aes(xintercept = meanRatio), linetype='dashed') +
  geom_vline(data = thresh, aes(xintercept = upRatio), linetype='dashed', colour='red') +
  geom_vline(data = thresh, aes(xintercept = downRatio), linetype='dashed', colour='blue') +
  geom_vline(data = multi.thresh, aes(xintercept = meanRatio), linetype='dotted') + 
  geom_vline(data = multi.thresh, aes(xintercept = upRatio), linetype='dotted', colour='pink') +
  geom_vline(data = multi.thresh, aes(xintercept = downRatio), linetype='dotted', colour='lightblue') +
  labs(title = "Genotype consistency")

plot3 <- ggplot(by.var_library %>% filter(total < 1000), mapping=aes(x=total)) + geom_density() + 
  facet_wrap(vars(library, gt.label), scales='free') +
  geom_vline(data = thresh, aes(xintercept = meanCov), linetype='dashed') +
  geom_vline(data = thresh, aes(xintercept = upCov), linetype='dashed', colour='red') +
  geom_vline(data = thresh, aes(xintercept = downCov), linetype='dashed', colour='blue') +
  geom_vline(data = multi.cov, aes(xintercept = meanCov), linetype='dotted') + 
  geom_vline(data = multi.cov, aes(xintercept = upCov), linetype='dotted', colour='pink') +
  geom_vline(data = multi.cov, aes(xintercept = downCov), linetype='dotted', colour='lightblue') +
  labs(title = "Sum of iPCR counts for all fragments")

plot4 <- ggplot(by.var_library, mapping=aes(x=okprop)) + geom_density() + 
  facet_wrap(vars(library, gt.label), scales='free') +
  geom_vline(data = thresh, aes(xintercept = meanOk), linetype='dashed') +
  geom_vline(data = thresh, aes(xintercept = downOk), linetype='dashed', colour='blue') +
  labs(title = "Proportion of non-imputed fragments with allele expected from genotype")

pdf(file=paste0(out_pref, ".threshold_plots.multi.pdf"), width=11)
print(plot1)
print(plot2)
print(plot3)
print(plot4)
dev.off()

message(sprintf("Done plots, starting to do filtering: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Start marking filter pass/fails
by.var_library <- by.var_library %>% left_join(thresh)
# Allele proportions
by.var_library <- by.var_library %>% 
  mutate(FltPropref=if_else(refPropEquiv >= downRatio & refPropEquiv <= upRatio, "PASS","FAIL")) %>% 
  replace_na(list(FltPropref='FAIL'))
# Total input
by.var_library <- by.var_library %>% 
  mutate(FltInput = if_else(total < 1000 & total >= downCov & total <= upCov, "PASS", "FAIL"))
# Imputed proportion
by.var_library <- by.var_library %>%
  mutate(FltImp = if_else(maxImp <= upImp, "PASS", "FAIL"))
# Remove groups if the proportion of 'other' is too high (based on mean and SD of ok prop)
by.var_library <- by.var_library %>%
  mutate(FltOK = if_else(okprop >= downOk, "PASS","FAIL"))

# Get fragments from only those that pass the filters
filter_details <- by.var_library %>% select(-starts_with(c("up","down","mean")))
by.var_library <- by.var_library %>% filter(FltPropref == "PASS", FltInput == "PASS", FltImp == "PASS", FltOK == "PASS")

# Then as for biallelics remove alleles that aren't in the genotypes
counts.flt_all <- inner_join(counts, by.var_library %>% dplyr::select(VarID, library)) %>% filter(which.allele!='other')

message(sprintf("Done variant filtering, doing fragment filtering: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Filter the fragments list based on this
flt.frags.initial <- multi_data %>% inner_join(counts.flt_all %>% select(VarID:genotype))

## Add additional per-fragment filters here as required -- MAPQ, imputed etc
if (minMAPQ) {
  flt.frags.initial <- flt.frags.initial %>%
    filter(vapply(strsplit(MAPQs, ",", fixed = TRUE),
                  function(x) min(as.integer(x)), 
                  integer(1),
                  USE.NAMES = FALSE) > minMAPQ)
}
if (filter.cav) {
  flt.frags.initial <- flt.frags.initial %>% filter(!(VarCav %in% cav.excl))
}

nfrags_init <- multi_data %>% group_by(VarID, library) %>% summarise(fragcount_orig=n(), .groups = 'drop')
filtered <- flt.frags.initial %>% group_by(VarID, library) %>% summarise(fragcount_filt=n(), .groups = 'drop')
filter_details <- filter_details %>% left_join(nfrags_init) %>% left_join(filtered)

message(sprintf("Done most filtering, checking allele counts: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Check that there is more than one allele for the variant across both libraries after all filtering
countsC <- flt.frags.initial %>% 
  pivot_wider(id_cols = VarID, names_from = library, values_from = genotype, values_fn = dplyr::first)
varsToKeep <- countsC %>% mutate(across(everything(), ~ str_remove_all(., "/"))) %>% 
  pivot_longer(cols = -VarID, values_drop_na = TRUE) %>% 
  separate_rows(value, sep = "") %>% filter(value != "") %>% select(-name) %>% 
  distinct() %>% dplyr::count(VarID) %>% filter(n > 1) %>% pull(VarID)

## Split into variants that are now (or were originally!) biallelic, and variants with still >2 alleles
flt.frags <- flt.frags.initial %>% filter (VarID %in% varsToKeep)
alleles <- flt.frags %>% ungroup %>% group_by(VarID) %>% 
  summarise(n_alleles=n_distinct(VarAllele))
varsToKeep.multi <- alleles %>% filter (n_alleles > 2) %>% pull(VarID)
varsToKeep.biall <- setdiff(varsToKeep, varsToKeep.multi)

message(sprintf("Finalising filtering details and outputting: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Work out which variants were excluded at this stage and output the filter details
incountsC <- countsC %>% pull(VarID) %>% unique()
failed <- setdiff(incountsC, varsToKeep)
# Convert to data.table if not already
setDT(filter_details)
# Create lookup tables for faster joins
still_multi_dt <- data.table(VarID = varsToKeep.multi, status = "PASS")
biall_dt <- data.table(VarID = varsToKeep.biall, status = "2ALL")
incounts_dt <- data.table(VarID = failed, status = "FAIL")
# Use joins instead of %in%
filter_details[, FltMultiAll := NA_character_]
# Join in order of precedence (most specific first)
filter_details[still_multi_dt, FltMultiAll := "PASS", on = "VarID"]
filter_details[biall_dt, FltMultiAll := "2ALL", on = "VarID"]
filter_details[incounts_dt, FltMultiAll := "FAIL", on = "VarID"]

fwrite(filter_details, file=paste0(out_pref, ".filter_details_multi.tsv"), sep = "\t", 
quote = FALSE, row.names = FALSE)

message(sprintf("After filtering to vars to keep: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Process into chunks for model testing, using data.table
# Keep flt.frags as tibble to help with last section below
flt.frags.final <- flt.frags %>% filter(VarID %in% varsToKeep.multi)
setDT(flt.frags.final)
# Encode allele as factor
flt.frags.final[ , VarAllele := as.factor(VarAllele)]

output_frags <- function(frags) {
  # Make into long format
  dat <- melt(
    frags,
    id.vars = c("VarID","VarAllele", "input", "library","Barcode"), # keep these
    measure.vars = reps_rename, # columns to pivot
    variable.name = "cells", # new column holding the old name
    value.name = "cDNA", # new column holding the numeric count
    variable.factor = FALSE # keep 'cells' as plain character
  )

  # split "cells" into two new columns: species (Cow / Hum) and Rep (Rep1…)
  dat[, c("cellSpecies", "Rep") := tstrsplit(cells, "_", fixed = TRUE, keep = c(1L, 2L))]
  dat[, cells := NULL]

  return(dat)
}

# Get chunks and write each to file
dat_multi <- output_frags(flt.frags.final)

# Get unique VarIDs and assign to batches
# Use twice as many chunks as for main data set
chunk_count_multi <- 2*chunk_count
unique_var_ids <- unique(dat_multi$VarID)
batch_size <- ceiling(length(unique_var_ids) / chunk_count_multi)
var_id_batches <- cut(seq_along(unique_var_ids), 
                      breaks = chunk_count_multi, 
                      labels = FALSE)
batch_lookup <- data.table(VarID = unique_var_ids, batch = var_id_batches)

# Write each batch without modifying the original data
for(batch_num in 1:chunk_count_multi) {
  batch_var_ids <- batch_lookup[batch == batch_num, VarID]
  batch_data <- dat_multi[VarID %in% batch_var_ids]
  fwrite(batch_data, 
         file = paste0("chunk_multi.", batch_num, ".tsv.gz"),
         sep = "\t")
  
  message(sprintf("Written batch %d/%d (%d rows, %d unique VarIDs)", 
                  batch_num, chunk_count_multi, nrow(batch_data),
                  length(unique(batch_data$VarID))))
  rm(batch_data)
  gc()
}

message(sprintf("Done multi data, now processing biallelics: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

rm(flt.frags.final)
gc()

# Output the variants that are now biallelic and can be tested in the main model script,
# recoding as 0 and 1
# Output as additional chunk for each chromosome
# Filter the fragments down to the final set and recode alleles as 0 and 1
flt.frags.biall <- flt.frags %>% filter(VarID %in% varsToKeep.biall)
rm(flt.frags)
gc()
biall.byvar <- flt.frags.biall %>% group_by(VarID, VarAllele) %>% 
  group_keys() %>% ungroup() %>% group_by(VarID) %>% arrange(VarID, VarAllele) %>%
  mutate(NewAllele = row_number()-1)

flt.frags.biall <- flt.frags.biall %>% left_join(biall.byvar) %>% 
  select(-VarAllele) %>% rename(VarAllele='NewAllele') %>% 
  relocate(VarAllele, .after=VarID)

setDT(flt.frags.biall)
# fails if none on a particular chromosome
get_onechr <- function(chromosome) {
  chrfrags <- flt.frags.biall[.(chromosome), on = .(chr), nomatch = 0L]
  chr_output <- output_frags(chrfrags)
}
chr_newchunks <- lapply(chr, get_onechr)
for(j in seq_along(chr_newchunks)) {
  fwrite(chr_newchunks[[j]], 
         file = paste0("chunk.",chr[j], "_", chunk_count+1, ".tsv.gz"),
         sep = "\t")
}

# Save the recoding information
write_tsv(biall.byvar, file=paste0(out_pref, ".multi_now_bi.recoding.tsv"))


