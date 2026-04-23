#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))

# Get the library sizes and calculate library-specific normalisation factors to multiply by for each replicate
# Need to do each library separately as columns may not be in the same order
# Currently calculating normalisation factor as iPCR/cDNA so that when multiplied by the fragment
# values we are doing cDNA/iPCR for the fragment divided byt cDNA/iPCR for the library
# (different from previous version)

# All library size files should be in working directory

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
out_pref <- args[2]

source(config_file)

libsizes_list <- list()
for (i in (1:length(libs))) {
  libsize_files <- list.files(".", pattern = paste0("^", libs[i], ".*\\.lib_counts.txt"))
  libsizes1 <- read_tsv(libsize_files, show_col_types = FALSE)
  libsizes1 <- libsizes1 %>% mutate(library= libs[i])
  libsizes1 <- libsizes1 %>% group_by(library) %>% summarise_all(sum)
  libsizes_list[[i]] <- libsizes1
}
libsizes <- bind_rows(libsizes_list)

# Rename count columns to standardised format
libsizes <- libsizes %>% rename(!!!setNames(reps,reps_rename))

# Calculate normalisation factors using same iPCR for every Rep
cdnas_norm <- libsizes %>% mutate (across(all_of(reps_rename), ~ iPCR / .), .keep = 'unused') %>%
 relocate (library, .after = last_col())

# Make into long format and split Rep names to match chunks file
cdnas_norm <- cdnas_norm %>% 
  pivot_longer(cols = -library, names_to = "Experiment", values_to = "normf") %>%
  separate_wider_delim(Experiment, "_", names = c('cellSpecies','Rep'))

write_tsv(cdnas_norm, file=paste0(out_pref, ".cDNA_normfactors.tsv"))
