#!/usr/bin/env Rscript

# Calculate mean SuRE scores and skewness for all multiallelic variants
# Only uses fragments and variant/animal combinations that went into the model testing
# Will compute a score for each different cellSpecies and one for all data

suppressMessages(library(data.table))
suppressMessages(library(e1071))

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
normfactors <- args[2]
out_pref <- args[3]
source(config_file)

# Set threads for data.table as a parameter given to the script
n_threads <- as.integer(args[4])
setDTthreads(n_threads)

# Load the normalisation factors
cdnas_norm <- data.table::fread(normfactors, sep = "\t")
data.table::setkey(cdnas_norm, library, cellSpecies, Rep)

# Find and load the chunks files in this directory
chunkfiles <- list.files(pattern = "chunk_multi")
chunklist <- lapply(chunkfiles, fread)
chunks_all <- data.table::rbindlist(chunklist)
chunks_all[,Barcode := NULL]
rm(chunklist)
gc()

# Normalise the cDNA counts by library and transfection
chrdata <- cdnas_norm[chunks_all, on =.(library, cellSpecies, Rep)]       # left join chunks_all with norms
rm(chunks_all)
gc()
chrdata[, SuREScore := (cDNA / input) * normf] # compute score in-place

# Fragment counts per allele
by_var_frags <- chrdata[, .N, by = .(VarID, VarAllele)]
by_var_frags <- dcast(by_var_frags, VarID ~ VarAllele, value.var = "N", fill = 0)
setnames(by_var_frags, setdiff(names(by_var_frags), "VarID"), paste0("count_", setdiff(names(by_var_frags), "VarID")))

# Mean score by allele and species (cow/hum/all)
means_spec <- chrdata[, .(mean = mean(SuREScore)), by = .(spec=cellSpecies, VarID, VarAllele)]
means_all <- chrdata[, .(spec="all", mean = mean(SuREScore)), by = .(VarID, VarAllele)]
means_comb  <- rbindlist(list(means_spec, means_all), use.names = TRUE)

means_wide <- dcast(
  means_comb,
  VarID ~ spec + VarAllele,
  value.var = "mean",
  fill = NA_real_
)

rm(means_spec, means_all, means_comb)
gc()

# Skew
skew_spec <- chrdata[, .(skew = e1071::skewness(SuREScore, na.rm = TRUE, type = 3)),
                   by = .(spec = cellSpecies, VarID)]
skew_all <- chrdata[, .(spec = "all", skew = e1071::skewness(SuREScore, na.rm = TRUE, type = 3)),
                   by = .(VarID)]
skew_comb <- rbindlist(list(skew_spec, skew_all), use.names = TRUE)
skew_wide  <- dcast(skew_comb, VarID ~ spec, value.var = "skew", fill = NA_real_)
setnames(skew_wide, setdiff(names(skew_wide), "VarID"), paste0(setdiff(names(skew_wide), "VarID"), "_skew"))

rm(skew_spec, skew_all, skew_comb)
gc()

# Assemble final output
out <- merge(by_var_frags, means_wide, by = "VarID", all.x = TRUE)
out <- merge(out, skew_wide, by = "VarID", all.x = TRUE)

# Write output
fwrite(out, file = paste("multi", out_pref, "meanScores.tsv.gz", sep="."), sep = "\t")
