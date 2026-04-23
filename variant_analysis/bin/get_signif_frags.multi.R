#!/usr/bin/env Rscript

start_time <- Sys.time()

suppressMessages(library(data.table))

# Same as get_signif_frags but for (true) multiallelic variants
# Currently inefficient as loading all chunks for every chr, but getting away with
# it as chunks are small

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

resfile <- paste0(out_pref, ".multiallelic.tsv.gz")
message("got filename")
vars <- fread(resfile)
message("Loaded vars")
vars <- vars[Chr == chrid]
message("Filter to this chr")
pval_cols <- grep("p\\.adjust$", names(vars), value = TRUE)
message("got pval_cols")
vars <- vars[rowSums(vars[, ..pval_cols] <= pval_cutoff, na.rm = TRUE) > 0]
message("extracted signif")
signif <- vars[,c("VarID")]
message("got signif IDs")

cdnas_norm <- data.table::fread(normfactors, sep = "\t")
message("loaded normfactors")
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
  if (nrow(signif_frags)) {
    frag_list_orig[[i]] <- data.table::copy(signif_frags)
    signif_frags <- merge(signif_frags, cdnas_norm, by = c("library", "cellSpecies", "Rep"))
    signif_frags[, SuREScore := (cDNA / input) * normf]
    signif_frags <- signif_frags[, -c("normf","input","cDNA")]
    
    signif_frags <- dcast(signif_frags, VarID + VarAllele + Barcode + library  ~ cellSpecies+Rep,
                          value.var = c("SuREScore")
    )
    frag_list[[i]] <- signif_frags
  }  else {
    message(paste0("No fragments from chunk ", i))
  }
}

message(sprintf("Found the fragments: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

frags <- rbindlist(frag_list)
# Output fragment list
if (nrow(frags) == 0) {
  message (paste0("No fragments found for chr ", chrid))
  blank <- "No_fragments"
  write(blank, file = paste(chrid, out_pref, "multi.signif_frag_scores.txt.gz", sep="."), sep = "\t")
  write(blank, file = paste(chrid, out_pref, "multi.signif_frag_chunks.txt.gz", sep="."), sep = "\t")
  write(blank, file = paste(chrid, out_pref, "multi.signif_frag_positions.txt.gz", sep="."), sep = "\t")
  quit()
}

fwrite(frags, file = paste(chrid, out_pref, "multi.signif_frag_scores.txt.gz", sep="."), sep = "\t")

frags_orig <- rbindlist(frag_list_orig)
fwrite(frags_orig, file = paste(chrid, out_pref, "multi.signif_frag_chunks.txt.gz", sep="."), sep = "\t")
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

fwrite(positions, file = paste(chrid, out_pref, "multi.signif_frag_positions.txt.gz", sep="."), sep = "\t")

message(sprintf("Finished: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
