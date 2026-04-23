#!/usr/bin/env Rscript

start_time <- Sys.time()

suppressMessages(library(data.table))

# Extract fragments for significant variants, for plotting
# Finds the barcodes of the relevant fragments, normalises fragment scores, and produces two outputs:
# Table with the data needed for plotting except fragment position - VarAllele, library, normalised scores
# Get variant-level data (skew, p-value) during plotting step, not saved here
# List of barcodes and their positions, for processing at the next step to get x-axis positions

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

# process extracts snps for just this chr but doesn't bother for non-SNPs
snpsfile <- paste0(chrid,".snps.tsv")
otherfile <- paste0(out_pref,".non_snps.tsv.gz")

snps <- fread(snpsfile)
snps <- snps[Chr == chrid]
pval_cols <- grep("p\\.adjust$", names(snps), value = TRUE)
snps <- snps[rowSums(snps[, ..pval_cols] <= pval_cutoff, na.rm = TRUE) > 0]
signif <- snps[,c("VarID")]

message(sprintf("Loaded SNPs: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

other <- fread(otherfile)
other <- other[Chr == chrid]
other <- other[rowSums(other[, ..pval_cols] <= pval_cutoff, na.rm = TRUE) > 0]
signif2 <- other[,c("VarID")]
signif <- rbind(signif, signif2)

message(sprintf("Loaded others: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# to add when multi data is ready
# multifile <- paste0(out_pref, ".multi.tsv.gz")
# multi <- fread(multifile)
# multi <- multi[Chr == chrid]
# multi <- multi[rowSums(multi[, ..pval_cols] <= pval_cutoff, na.rm = TRUE) > 0]
# signif3 <- multi[,c("VarID")]
# signif <- rbind(signif, signif3)

# Load the normalisation factors
cdnas_norm <- data.table::fread(normfactors, sep = "\t")
data.table::setkey(cdnas_norm, library, cellSpecies, Rep)

message(sprintf("Starting to go through chunks files: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Get fragments from the chunks and calculate their normalised scores
frag_list <- list()
frag_list_orig <- list() # for just a subset of the original chunks files
chunkfiles <- list.files(pattern = "chunk")
for (i in seq_along(chunkfiles)) {
  chunkdata <- fread(chunkfiles[i])
  message(paste0("Loaded chunk ", i))
  signif_frags <- merge(chunkdata, signif, by = "VarID")
  frag_list_orig[[i]] <- data.table::copy(signif_frags)
  signif_frags <- merge(signif_frags, cdnas_norm, by = c("library", "cellSpecies", "Rep"))
  signif_frags[, SuREScore := (cDNA / input) * normf]
  signif_frags <- signif_frags[, -c("normf","input","cDNA")]
  
  signif_frags <- dcast(signif_frags, VarID + VarAllele + Barcode + library  ~ cellSpecies+Rep,
                value.var = c("SuREScore")
                )
  frag_list[[i]] <- signif_frags
}

message(sprintf("Found the fragments: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

frags <- rbindlist(frag_list)
# Output fragment list
if (nrow(frags) == 0) {
  message ("No fragments found, problem")
  quit(status = 1)
}

fwrite(frags, file = paste(chrid, out_pref, "signif_frag_scores.txt.gz", sep="."), sep = "\t")

frags_orig <- rbindlist(frag_list_orig)
fwrite(frags_orig, file = paste(chrid, out_pref, "signif_frag_chunks.txt.gz", sep="."), sep = "\t")
rm(frags_orig)
rm(frag_list_orig)

# Use the barcodes to find these fragments in the var counts tables and output
barcodes <- frags[,c("Barcode","library")]
barcodes <- unique(barcodes)
rm(frags)

message(sprintf("Found the barcodes, starting to go through counts tables: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

keepcols <- c("Barcode", "Start","End","Haplotype")
pos_list <- list()
for (i in seq_along(libs)) {
  resfile <- paste0(libs[i], ".", chrid, ".var_counts_table.txt.gz")
  resdata <- fread(resfile)[, ..keepcols]
  message(paste0("Loaded library ", resfile))
  resdata <- unique(resdata)
  resdata <- merge(resdata, barcodes)
  pos_list[[i]] <- resdata
}
positions <- rbindlist(pos_list)

fwrite(positions, file = paste(chrid, out_pref, "signif_frag_positions.txt.gz", sep="."), sep = "\t")

message(sprintf("Finished: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))