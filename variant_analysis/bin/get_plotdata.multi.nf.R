#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(purrr))

# Combine test results from all multi chunks and for each transfected cell species
# Remove variants with warnings in either of the comparison models
# Correct pvalues
# Add variant details (skews, normalised expression, genotypes and variant types)

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

# Set threads for data.table as a parameter given to the script
n_threads <- as.integer(args[3])
setDTthreads(n_threads)

source(config_file)

drop_cols <- c("var_warnings","novar_warnings")
# Code for getting the relevant data
get_data <- function(spec, filestem) {
  files <- list.files (pattern = paste0(filestem,"[[:digit:]]+.", spec, ".tsv.gz"))
  data <- rbindlist(lapply(files, fread))
  
  # remove vars with warnings in either of the two models compared
  warned <- data[(var_warnings | novar_warnings)]
  data <- data[!(var_warnings | novar_warnings)]
  
  # Rename columns (except VarID)
  old_names <- setdiff(names(data), "VarID")
  new_names <- paste0(spec, "_", old_names)
  new_names[3] <- paste0(spec, "_", "p.value")
  data <- data[, -drop_cols, with = FALSE]
  setnames(data, old_names, new_names, skip_absent = TRUE)
  setnames(warned, old_names, new_names, skip_absent = TRUE)
  
  return(list(warned, data))
}

# Get the model results
multi_filestem <- "modelResults.chunk_multi."
  
dataList <- lapply(specs, get_data, filestem=multi_filestem)
warnList <- lapply(dataList, `[[`, 1)
warnings <- reduce(warnList, merge, by="VarID", all=TRUE)
fwrite(warnings, file=paste0(out_pref, '.multiallelic.model_warnings.tsv.gz'), sep="\t", na = "NA")

multiList <- lapply(dataList, `[[`, 2)
rm(dataList)
multi <- reduce(multiList, merge, by="VarID", all=TRUE)
rm(multiList)

# Do p-value correction for each column
for (species in specs) {
  pval_orig <- paste0(species, "_p.value")
  pval_correct <- paste0(species, "_p.adjust")
  multi[, (pval_correct) := p.adjust(get(pval_orig), method = "BH")]
}

# Load vardetails to get chromosome and genotypes
vardataList <- list()
for (i in seq_along(chr)) {
  chrid <- chr[i]
  chr_varfile <- paste0(chrid, ".", out_pref, ".var_details.tsv")
  chr_vardata <- fread(chr_varfile)
  chrinfo <- merge(multi[,"VarID"], chr_vardata)
  vardataList[[i]] <- chrinfo
}
vardata <- rbindlist(vardataList)
multi <- merge(multi, vardata, by ="VarID", all.x = T)

# Add alleles
more_vardata <- fread(paste0(out_pref, ".", "alleles.tsv"), col.names = c("VarID", "Pos","RefBase","AltBase"))
multi <- merge(multi, more_vardata, by="VarID", all.x = T)

# add scores
scores <- fread(paste("multi", out_pref, "meanScores.tsv.gz", sep="."))
multi <- merge(multi, scores, by = "VarID", all.x = T)

check <- multi[is.na(RefBase) | is.na(all_skew) | is.na(Chr)]
if (nrow(check)) {
  message <- "Some variants were tested but don't have the expected details"
  stop(message)
}

# Write output
fwrite(multi, file=paste0(out_pref, ".multiallelic.tsv.gz"), sep="\t", na = "NA")
