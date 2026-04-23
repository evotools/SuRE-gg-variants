#!/usr/bin/env Rscript

# Calculate mean SuRE scores and skewness for each variant
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
chrid <- as.character(args[2])
normfactors <- args[3]
out_pref <- args[4]
source(config_file)

# Set threads for data.table as a parameter given to the script
n_threads <- as.integer(args[5])
setDTthreads(n_threads)

# Load the normalisation factors
cdnas_norm <- data.table::fread(normfactors, sep = "\t")
data.table::setkey(cdnas_norm, library, cellSpecies, Rep)

# Find and load the chunks files in this directory
chunkfiles <- list.files(pattern = "chunk")
chunklist <- lapply(chunkfiles, fread)
#data <- do.call(c, chunklist)
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
setnames(by_var_frags, c("0","1"), c("count_0","count_1"))

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

# Find which libraries are in the output and add genotypes for only libraries in the analysis
# For variants which were originally non-biallelic (but now are), add the original genotypes
libs_byvar <- unique(chrdata[, .(VarID, library)])
chr_varfile <- paste0(chrid, ".", out_pref, ".var_details.tsv")
chr_vardata <- fread(chr_varfile)
vardata_bylib <- melt(
  chr_vardata,
  id.vars = "VarID",
  measure.vars = libs,
  variable.name = 'library',
  value.name = 'genotype'
)
libs_byvar <- merge(libs_byvar, vardata_bylib, by = c('VarID', 'library'), all.x = T)
var_gt <- dcast(
  libs_byvar,
  VarID  ~ library,
  value.var = "genotype"
)

# Assemble final output
out <- merge(by_var_frags, means_wide, by = "VarID", all.x = TRUE)
out <- merge(out, skew_wide, by = "VarID", all.x = TRUE)
out <- merge(out, var_gt, by = "VarID", all.x = TRUE)
out[, Chr := chrid]
setcolorder(out, c("VarID","Chr"))

# Write output
fwrite(out, file = paste(chrid, out_pref, "meanScores.tsv.gz", sep="."), sep = "\t")
