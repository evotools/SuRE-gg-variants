#!/usr/bin/env Rscript

start_time <- Sys.time()

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# Take the variant-labelled fragments for each chromosome (combined in combine_libs.R), filter for various attributes and process into the data chunks to use for the model test
# Instead of keeping count of what's filtered out as you go along, mark whether the variants pass the filter or not

# Stop on any error immediately
options(error = function() {
  cat("Error occurred at:", Sys.time(), "\n")
  cat("Error message:", geterrmessage(), "\n")
  traceback()
  quit(status = 1)  # Exit with error code
})
options(warn = 2)  # Turn warnings into errors

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
chrid <- as.character(args[2])
resfile <- args[3]
out_pref <- args[4]
source(config_file)

# Set threads for data.table as a parameter given to the script
n_threads <- as.integer(args[5])
setDTthreads(n_threads)

message(sprintf("Before loading data: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
# Get fragment data for this chromosome
onechr <- readRDS(resfile)

message(sprintf("After loading data: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Get the possible combinations of VarID + library, and their associated genotypes, iPCR totals, refprops and impprops
# If no data on imputed fragments then add a dummy column of 0s
if ("impProp" %in% colnames(onechr)) {
  has.imp <- T
  onechr <- onechr %>% group_by(VarID, library, genotype) %>% mutate(maxImp=max(impProp), .keep = 'unused') %>% ungroup
} else{
  has.imp <- F
  onechr <- onechr %>% ungroup %>% mutate(maxImp=0)
}
by.var_library <- onechr %>% group_by(VarID, library, genotype, varInput, propRef, maxImp) %>% summarise(fragcount_orig=n(), .groups = 'drop')

# Load in thresholds calculated across all chr and add to variants object
ratioThresh <- read_tsv(paste0(out_pref,".thresholds.ratio.tsv"), show_col_types = F)
covThresh <- read_tsv(paste0(out_pref,".thresholds.iPCR.tsv"), show_col_types = F)
by.var_library <- by.var_library %>% left_join(ratioThresh) %>% left_join(covThresh)

# Add filter for variants where our genotypes dont seem to agree with the relative frequencies of fragments carrying the different alleles
by.var_library <- by.var_library %>% mutate(FltPropref=if_else(propRef >= downRatio & propRef <= upRatio, "PASS","FAIL")) %>% replace_na(list(FltPropref='FAIL'))

# Add filter for variants with unusual coverage with iPCR fragments. To get thresholds, in a previous step did a hard filter for >1000, then find +/- 2.58*SD from mean and exclude if outside this. Use the precalculated thresholds and re-apply the <1000 threshold (just in case)
by.var_library <- by.var_library %>% mutate(FltInput = if_else(varInput < 1000 & varInput >= downCov & varInput <= upCov, "PASS", "FAIL"))

# If imputation data available, filter out variants with too high proportion of imputed fragments
# For hets, filter out the variant+library if the maximum of the two imputation proportions (calculated earlier) is >cutoff
if (has.imp) {
  impThresh <- read_tsv(paste0(out_pref,".thresholds.imp.tsv"), show_col_types = F)
  by.var_library <- by.var_library %>% left_join(impThresh)
  by.var_library <- by.var_library %>% mutate(FltImp = if_else(maxImp <= upImp, "PASS", "FAIL"))
} else {
  by.var_library <- by.var_library %>% mutate(FltImp = PASS)
}

# Get fragments from those that pass the filters
filter_details <- by.var_library %>% select(-starts_with(c("up","down","mean")))
by.var_library <- by.var_library %>% filter(FltPropref == "PASS", FltInput == "PASS", FltImp == "PASS")

# Those that are left have fragment ratios consistent with our genotypes but due to sequencing errors etc we potentially still have fragments carrying the wrong allele at homozygote sites. So just drop these so don't confuse the models.
by.var_ani_all <- onechr %>% distinct(across(c(VarID, VarAllele, library, genotype)))
by.var_ani_all <- by.var_ani_all %>% inner_join(by.var_library %>% select(VarID, library, genotype))
by.var_ani_all <- by.var_ani_all %>% filter(!(genotype=="0/0" & VarAllele==1)) %>% filter(!(genotype=="1/1" & VarAllele==0))

message(sprintf("After filtering genotypes: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Filter the fragments list based on this
flt.frags <- onechr %>% select(-propRef, -varInput, -maxImp) %>% inner_join(by.var_ani_all)
rm(onechr)
gc()

## Add additional per-fragment filters here as required -- MAPQ, imputed etc
if (minMAPQ) {
  flt.frags <- flt.frags %>%
    filter(vapply(strsplit(MAPQs, ",", fixed = TRUE),
                  function(x) min(as.integer(x)), 
                  integer(1),
                  USE.NAMES = FALSE) > minMAPQ)
}
if (filter.cav) {
  flt.frags <- flt.frags %>% filter(!(VarCav %in% cav.excl))
}

filtered <- flt.frags %>% group_by(VarID, library, genotype) %>% summarise(fragcount_filt=n(), .groups = 'drop')
filter_details <- filter_details %>% left_join(filtered)

message(sprintf("After filtering fragments: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Filter to sites (based on GT) where we see both alleles across the two libraries (could be opposing homozygotes). Filtered for this at the start of the pipeline but this may have changed after filtering (e.g. need to remove variants where one library is hom and the other has been removed by the filters above)
# This allows a variant to stay if one library is 0/1 and the other is NA (i.e. has no fragments or has been filtered out) but needed to change the NA to another string (used 99). Just needs to be >1 allele in the whole dataset
# First calculate genotype table with filtered fragments
countsC <- flt.frags %>% pivot_wider(id_cols = VarID, names_from = library, values_from = genotype, values_fn = first)
countsC <- countsC %>% replace_na(list("BOAEC-D1"="99", "Bi-83-A"="99"))
varsToKeep<-countsC[which((str_count(countsC$`Bi-83-A`, "0")+str_count(countsC$`BOAEC-D1`, "0"))>0 & (str_count(countsC$`Bi-83-A`, "1")+str_count(countsC$`BOAEC-D1`, "1"))>0),] %>% drop_na() %>% pull(VarID)

# Work out which variants were excluded at this stage and output the filter details
incountsC <- countsC %>% pull(VarID) %>% unique()
failed <- setdiff(incountsC, varsToKeep)
setDT(filter_details)
passed_dt <- data.table(VarID = varsToKeep, status = "PASS")
failed_dt <- data.table(VarID = failed, status = "FAIL")
filter_details[, FltMultiAll := NA_character_]
# Join in order of precedence (most specific first)
filter_details[passed_dt, FltMultiAll := "PASS", on = "VarID"]
filter_details[failed_dt, FltMultiAll := "FAIL", on = "VarID"]
write_tsv(filter_details, file=paste(chrid, out_pref, "filter_details.tsv", sep="."))

message(sprintf("After filtering to vars to keep: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Process into chunks for model testing, using data.table
setDT(flt.frags)
setkey(flt.frags, VarID)
# Get final set
flt.frags.final <- flt.frags[.(varsToKeep), nomatch = 0L]
rm(flt.frags)
gc()

# Make into long format
dat <- melt(
  flt.frags.final,
  id.vars = c("VarID","VarAllele", "input", "library","Barcode"), # keep these
  measure.vars = reps_rename, # columns to pivot
  variable.name = "cells", # new column holding the old name
  value.name = "cDNA", # new column holding the numeric count
  variable.factor = FALSE # keep 'cells' as plain character
)
rm(flt.frags.final)
gc()

message(sprintf("Converted to long format: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# split "cells" into two new columns: species (Cow / Hum) and Rep (Rep1…)
dat[, c("cellSpecies", "Rep") := tstrsplit(cells, "_", fixed = TRUE, keep = c(1L, 2L))]
dat[, cells := NULL]

message(sprintf("Done split cellSpecies column, starting to write chunks: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Output this data table in batches -- now doing grouping just before running the 
# model to save time on I/O

# Get unique VarIDs and assign to batches
unique_var_ids <- unique(dat$VarID)
batch_size <- ceiling(length(unique_var_ids) / chunk_count)

# Create batch assignments
var_id_batches <- cut(seq_along(unique_var_ids), 
                      breaks = chunk_count, 
                      labels = FALSE)

# Create lookup table
batch_lookup <- data.table(VarID = unique_var_ids, batch = var_id_batches)

# Write each batch without modifying the original data
for(batch_num in 1:chunk_count) {
  batch_var_ids <- batch_lookup[batch == batch_num, VarID]
  
  # Filter data for this batch
  batch_data <- dat[VarID %in% batch_var_ids]
  
  fwrite(batch_data, 
         file = paste0("chunk.", chrid, "_", batch_num, ".tsv.gz"), 
         sep = "\t")
  
  message(sprintf("Written batch %d/%d (%d rows, %d unique VarIDs)", 
                  batch_num, chunk_count, nrow(batch_data),
                  length(unique(batch_data$VarID))))
  rm(batch_data)
  gc()
}

message(sprintf("After writing data (end of script): %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
