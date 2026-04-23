#!/usr/bin/env Rscript

start_time <- Sys.time()

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# Originally a script to normalise SuRE scores (so based on previous normalise_SuREdata.R scripts) but have now take this out
# Now combines fragments from multiple libraries without normalising, and removes fragments where the barcode appears in more than one library (likely misassignment during Illumina demultiplexing?)
# Then gets stats for number of variants, type, number of imputed, errors, etc

# Adapted to run within Nextflow pipeline

# Stop on any error immediately
options(error = function() {
  cat("ERROR at", Sys.time(), ":\n")
  cat("Message:", geterrmessage(), "\n")
  cat("Last call:", deparse(sys.calls()[[length(sys.calls())]]), "\n")
  cat("Traceback:\n")
  traceback(max.lines = 10, full = TRUE)
  quit(status = 1)
})
options(warn = 2)  # Turn warnings into errors

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
chrid <- as.character(args[2])
out_pref <- args[3]
gts_loc <- if (length(args) >= 4) args[4] else NA

source(config_file)
n_threads <- as.integer(args[5])
setDTthreads(n_threads)

message(sprintf("Before loading data: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Load in the data for this chromosome
# Rename the cDNA columns to standardise for later
onechr <- NULL
barcodes <- NULL
keepcols <- c("VarID", "Barcode", "VarAllele", "iPCR", reps, "VarType", "VarCav", "MAPQs")
for (i in seq_along(libs)) {
  resfile <- paste0(libs[i], ".", chrid, ".var_counts_table.txt.gz")
  resdata <- fread(resfile, colClasses = list(character = "MAPQs"))[, ..keepcols] %>% as_tibble()
  resdata <- resdata %>% rename(!!!setNames(c(reps,"iPCR"),c(reps_rename,"input")))
  resdata <- resdata %>% add_column(library=libs[i])
  bc_onelib <- resdata %>% select(Barcode, input, library) %>% distinct()
  barcodes <- bind_rows(barcodes, bc_onelib)
  onechr <- bind_rows(onechr,resdata)
}

message(sprintf("Loaded data: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# If analysing more than one library, find duplicated barcodes (appears in >1 library) and exclude
if (length(libs) > 1) {
  dup_bc <- barcodes %>% group_by(Barcode) %>% summarise(in_samples=n_distinct(library)) %>% filter (in_samples > 1)
  dup_bc_id <- dup_bc %>% pull(Barcode)
  onechr <- onechr %>% filter(!(Barcode %in% dup_bc_id))
}

message(sprintf("Removed duplicate barcodes: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Separate into biallelic and non-biallelic variants and output both for further processing.
onechr_multi <- onechr %>% filter(VarType == "multi")
saveRDS(onechr_multi,file=paste(chrid, out_pref, "multi.comb.rds", sep="."))

vartypes <- onechr %>% distinct(VarID, VarType)
write_tsv(vartypes, file=paste(chrid, out_pref, "vartypes.tsv", sep="."))

onechr <- onechr %>% filter(VarType != "multi") %>% select(-VarType)

message(sprintf("Removed multiallelics and output variant types: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Load in genotypes if appropriate (if GT file not given then assume all are het i.e. human library)
# GTs from cactus graph genome are phased, need to make sure 0|1 is equivalent to 1|0
if (!is.na(gts_loc))  {
  gts <- read_tsv(gts_loc, col_types = "cciccc") %>% select(-Chr, -RefLength, -AltLengths) %>% pivot_longer(cols=-VarID, names_to = "library", values_to = "genotype")
  onechr <- onechr %>% left_join(gts, join_by(VarID, library))
} else  {
  onechr <- onechr %>% mutate(genotype = "0/1")
}
# remove rows without genotypes (for if using a subset)
onechr <- onechr %>% drop_na(genotype)

message(sprintf("Loaded genotypes, starting to calculate filtering statistics: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Group into variants and calculate for each variant+library combination the proportion of reads that are ref, proportion that are imputed, and the total iPCR count for both libraries, for later filtering
# To check fragment ratios use only allele assignments obtained from the sequencing (i.e. not imputed). For total input use all
counts <- onechr %>% group_by(VarID, VarAllele, library, genotype) %>% summarise(totalInput=sum(input))

fs_counts <- onechr %>% filter (VarCav == 0 | VarCav == "E") %>% 
  group_by(VarID, VarAllele, library, genotype) %>% summarise(nonI=sum(input))
fragRatios <- fs_counts %>% full_join(counts) %>% 
  pivot_wider(names_from = VarAllele, values_from = c(nonI, totalInput)) %>% 
  replace_na(list(nonI_0 = 0, nonI_1 = 0, totalInput_0 = 0, totalInput_1 = 0)) %>%
  mutate(varnonI=nonI_0+nonI_1, propRef=nonI_0/varnonI, varInput=totalInput_0+totalInput_1, .keep = 'unused') %>% 
  select(-varnonI)
imputeCounts <- onechr %>% filter (!(VarCav == 0 | VarCav == "E")) %>% 
  group_by(VarID, VarAllele, library, genotype) %>% 
  summarise(totalImp=sum(input))
imputeProp_allele <- counts %>% left_join(imputeCounts) %>% 
  replace_na(list(totalImp = 0)) %>% 
  mutate(impProp = totalImp/totalInput, .keep='unused')
imputeProp_var <- counts %>% left_join(imputeCounts) %>% 
  replace_na(list(totalImp = 0)) %>% ungroup %>% 
  group_by(VarID, library, genotype) %>% 
  summarise(impProp = sum(totalImp)/sum(totalInput))
fragRatios_out <- fragRatios %>% left_join(imputeProp_var)

# Output fragRatios variant+library+genotype
saveRDS (fragRatios_out, file=paste(chrid, out_pref, "counts.rds", sep="."))

# Add them to fragment details (frags from same variant and allele have same value) and save that
onechr <- onechr %>% left_join(fragRatios) %>% left_join (imputeProp_allele)
saveRDS (onechr,file=paste(chrid, out_pref, "comb.rds", sep="."))
