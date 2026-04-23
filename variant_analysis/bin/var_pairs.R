#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(pscl))
suppressMessages(library(nplyr))
suppressMessages(library(parallel))
suppressMessages(library(furrr))
suppressMessages(library(broom))
suppressMessages(library(broom.helpers))
suppressMessages(library(lmtest))

# Go through pairs of variants that are found on the same fragments, and
# where there's a reasonable number of fragments and the two variants are not 
# in near-complete LD, test whether they have separate effects
# Find which of the two has the lower p-value and fit a model with the allele of
# this variant as a covariate, then compare with a model including both variants
# to see if this is a better fit
# Output all pairs of variants found on the same fragment, even if not tested (which will be most)

# Stop on any error immediately
options(error = function() {
  cat("Error occurred at:", Sys.time(), "\n")
  cat("Error message:", geterrmessage(), "\n")
  traceback()
  quit(status = 1)  # Exit with error code
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
chrid <- as.character(args[2])
out_pref <- args[3]
source(config_file)

n_threads <- as.integer(args[4])
setDTthreads(n_threads)

spec <- 'Cow' # assuming this is all we want, can set as input if needed

# lazy shortcut to avoid having to remove X earlier
if (!(chrid %in% chr)) {
  message (paste0("Skipping chr ", chrid))
  blank <- "No_fragments"
  write(blank, file = paste(out_pref, chrid, "variant_pairs.tsv.gz", sep="."))
  quit()
}

min_frags <- 10 # only test where there are at least 3 haplotypes with this many fragments

# Get variant results to find p-values for each significant (and exclude any where focal species wasn't signif)
snpsfile <- paste0(out_pref,".snps.tsv.gz")
otherfile <- paste0(out_pref,".non_snps.tsv.gz")
pval_col <- paste0(spec, "_p.adjust")
skewcol <- paste0(spec, "_skew")

snps <- fread(snpsfile)
snps <- snps[Chr == chrid]
snps <- snps[get(pval_col) <= pval_cutoff & get(skewcol) <= 20]
getcols <- c("VarID", pval_col)
snp.signif <- snps[, ..getcols]
rm(snps)

other <- fread(otherfile)
other <- other[Chr == chrid]
other <- other[get(pval_col) <= pval_cutoff & get(skewcol) <= 20]
other.signif <- other[,..getcols]
pvals <- rbind(snp.signif, other.signif)
rm(other)

# Load subset of chunks file containing only significant variants
sig_frags_long <- fread(paste0(chrid, ".all_cattle.signif_frag_chunks.txt.gz"))
sig_frags_long <- sig_frags_long[cellSpecies == spec]
sig_frags_long <- sig_frags_long[pvals[, .(VarID)], on = "VarID"]

# Find all the possible combinations on the same fragment
sig_frags <- unique(sig_frags_long[ ,.(VarID, Barcode, VarAllele, library)])
mergecols <- setdiff(colnames(sig_frags), c("VarID","VarAllele"))
sig_paired <- merge(sig_frags, sig_frags, by = mergecols, allow.cartesian = T, suffixes = c("_1","_2"))[VarID_1 != VarID_2]
sig_pairs <- unique(sig_paired[, .(VarID_1, VarID_2)])

# For each unique pair, find which has the smaller p-value -- this will be Var1 (lead variant)
# If a tie then test both
sig_pairs <- merge(sig_pairs, pvals, by.x = "VarID_1", by.y = "VarID", all.x = T)
sig_pairs <- merge(sig_pairs, pvals, by.x = "VarID_2", by.y = "VarID", suffixes = c("_1","_2"), all.x = T)
pval1 <- paste0(pval_col, "_1")
pval2 <- paste0(pval_col, "_2")
sig_pairs <- sig_pairs[ , keep := 
                          fifelse( get(pval1) < get(pval2), T,
                                   fifelse( get(pval1) == get(pval2), T, F)
                                   )]
testpairs <- sig_pairs[(keep), .(VarID_1,VarID_2)]

# Restrict to these (ordered) pairs
sig_paired <- merge(sig_paired, testpairs, by = c("VarID_1", "VarID_2"))

# Find counts and allele combinations
sig_paired <- sig_paired[, hap := paste0(VarAllele_1, VarAllele_2)]
pair_haps <- sig_paired[, .N, by = .(VarID_1, VarID_2, hap)]
pair_hap_counts <- dcast(pair_haps, VarID_1 + VarID_2 ~ hap, value.var = "N")
hapcols <- setdiff(colnames(pair_hap_counts), c("VarID_1", "VarID_2"))
pair_hap_counts <- pair_hap_counts[, haps := rowSums(!is.na(.SD) & .SD >= min_frags), .SDcols = hapcols]

modelpairs <- pair_hap_counts[haps >= 3]
modelpairs <- modelpairs[,.(VarID_1,VarID_2)]
modelpairs <- modelpairs[, pairID := paste(VarID_1, VarID_2, sep="_")]

# Get the fragments for just these pairs for the models and group by pair of vars
testvars <- data.table(VarID = unique(c(modelpairs$VarID_1, modelpairs$VarID_2)))
test_unpaired <- sig_frags_long[testvars, on = "VarID"]
#rm(sig_frags_long)
mergecols2 <- setdiff(colnames(test_unpaired), c("VarID","VarAllele"))
test_paired <- merge(test_unpaired, test_unpaired, by = mergecols2, allow.cartesian = T, suffixes = c("_1","_2"))
testfrags <- merge(test_paired, modelpairs, by = c("VarID_1","VarID_2"))
rm(test_paired)
testfrags <- testfrags[, pairID := paste(VarID_1, VarID_2, sep="_")]
testfrags <- testfrags[,-c("library", "Barcode", "cellSpecies", "VarID_1", "VarID_2")]

# Get into same format as for main model script
var_data <- split(testfrags, by = "pairID", keep.by = FALSE, flatten = FALSE)
#rm(by_var_flat)
#gc()
var_data <- tibble::tibble(pairID = names(var_data), data  = purrr::map(unname(var_data), tibble::as_tibble))

# Set up and run the models
modelv1 <- cDNA ~ offset(log(input)) + Rep + VarAllele_1 | 1
modelb <- cDNA ~ offset(log(input)) + Rep + VarAllele_1 + VarAllele_2 | 1

run_safe_quiet <- safely(quietly(function(whichmodel, df) {
    zeroinfl(whichmodel, data = df, dist = "negbin")
}))

modelprob <- function(safe.model) {
  failed <- !is.null(safe.model$error)
  warnings <- safe.model$result$warnings
  has.warning <- length(warnings) > 0
  return(failed || has.warning)
}

run_models <- function(df) {
  results <- list(onevar=NA, imp=NA, onevar_fail=F, both_fail=NA)
  
  var1result <- run_safe_quiet(modelv1, df)
  var1prob <- modelprob(var1result)
  
  if (var1prob) {
    results$onevar_fail <- T
    return(results)
  }
  
  # Get p-values for the one-variant model (might be subset of fragments)
  var1_pval <- tidy_zeroinfl(var1result$result$result) %>% filter(term == "VarAllele_1") %>% pull(p.value)
  results$onevar <- var1_pval
  
  results$onevar
  bothresult <- run_safe_quiet(modelb, df)
  bothprob <- modelprob(bothresult)
  if (bothprob) {
    results$both_fail <- T
    return(results)
  }
  results$both_fail <- F
  
  # Do both comparisons
  var2imp <- lrtest(bothresult$result$result, var1result$result$result)
  var2imp <- var2imp$`Pr(>Chisq)`
  var2imp <- var2imp[!is.na(var2imp)]
  results$imp <- var2imp
  
  return(results)
}

data_list <- var_data$data
comp_models <- mclapply(data_list, run_models, mc.cores = n_threads)
var_data$modelres <- comp_models
var_data <- var_data %>% select(-data)

var_data <- var_data %>% mutate(newVar1_pval= map_dbl(modelres, "onevar"),
                                twoall_better = map_dbl(modelres, "imp"),
                                failedv = map_lgl(modelres, "onevar_fail"),
                                failedb = map_lgl(modelres, "both_fail"), .keep = 'unused')

# Add to the list of all pairs of variants on the same fragment in either direction and whether this ordered combination was used
# reorder sig_pairs <- sig
sig_pairs <- merge(sig_pairs, pair_hap_counts, by = c("VarID_1","VarID_2"), all.x = T)
sig_pairs <- sig_pairs[, pairID := paste(VarID_1, VarID_2,sep="_")]
sig_pairs <- merge(sig_pairs, var_data, by=c("pairID"), all.x=T)

fwrite(sig_pairs, file=paste(out_pref, chrid, "variant_pairs.tsv.gz", sep="."), na = "NA", sep="\t", quote=F)


