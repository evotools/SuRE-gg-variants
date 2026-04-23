#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# For multiallelics that only had two alleles segregating and so were processed as biallelics,
# process the recoding and filtering data to get the alleles that were actually compared
# and label them with variant types based on these

# Stop on any error immediately
options(error = function() {
  cat("Error occurred at:", Sys.time(), "\n")
  cat("Error message:", geterrmessage(), "\n")
  traceback()
  quit(status = 1)  # Exit with error code
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
out_pref <- args[2]
source(config_file)

n_threads <- as.integer(args[3])
setDTthreads(n_threads)

alleles <- fread(paste0(out_pref, ".", "alleles.tsv"), col.names = c("VarID", "Pos","Ref","Alts"))
recoded <- fread(paste0(out_pref, ".multi_now_bi.recoding.tsv"))

recoded[,VarAllele := as.factor(VarAllele)]

# Get each allele in a separate column then join to recoded
multi_alleles <- alleles[VarID %in% recoded$VarID]
max_cols <- multi_alleles[, max(lengths(strsplit(Alts, ",", fixed = TRUE)))]
multi_alleles[, (as.character(1:max_cols)) := tstrsplit(Alts, ",", fixed = TRUE)]
setnames(multi_alleles, "Ref","0")

allele_cols <- as.character(0:max_cols)
alleles_long <- melt(
  multi_alleles,
  id.vars = "VarID",
  measure.vars = allele_cols,
  variable.name = 'VarAllele',
  value.name = 'Seq'
)
alleles_long <- na.omit(alleles_long)
both <- merge(recoded, alleles_long, by = c("VarID", "VarAllele"), all.x = T)

# Compare the two retained alleles
byvar <- dcast(
  both,
  VarID ~ VarAllele,
  value.var = "Seq"
)

#byvar <- merge(byvar, refs, by="VarID", all.x = T)

class_var <- function(alls_list) {
#  refall <- alls_list[[1]]
#  alls_list <- alls_list[-1]
  
  alls <- unlist(alls_list)
  kept <- alls[!is.na(alls)]
  
  # Find minimum length
  min_len <- min(nchar(kept[1]), nchar(kept[2]))
  if (is.na(min_len)) {
    message <- paste("Error:", kept[1], kept[2], sep = " ")
    return(message)
  }
  
  # Find common prefix length
  common_len <- 0
  for (i in 1:min_len) {
    if (substr(kept[1], i, i) == substr(kept[2], i, i)) {
      common_len <- i
    } else {
      break
    }
  }
  
  # Reverse to find common suffix length in the remaining sequence
  rev1 <- sapply(lapply(strsplit(kept[1], NULL), rev), paste, collapse = "")
  rev2 <- sapply(lapply(strsplit(kept[2], NULL), rev), paste, collapse = "")
  common_suff <- 0
  if (common_len < min_len) {
    for (i in 1:(min_len - common_len)) {
      if (substr(rev1, i, i) == substr(rev2, i, i)) {
        common_suff <- 0 - i
      } else {
        break
      }
    }
  }
  
  if ((common_suff - common_len) > min_len) {
    message <- paste("Problem with comparing alleles: ", kept[1], kept[2], common_len, common_suff)
    stop(message)
  }
  
  # For an indel in standard notation this would reduce one allele to "" so add this back
  if ((min_len - common_len + common_suff) < 1) {
    common_len <- common_len - 1
  }
  
  # Separate into shared prefix, reduced alleles and shared suffix
  prefix <- substr(kept[1], 1, common_len)
  suffix <- stringr::str_sub(kept[1], common_suff)
  
  all1 <- stringr::str_sub(kept[1], common_len+1, common_suff-1)
  all2 <- stringr::str_sub(kept[2], common_len+1, common_suff-1)
  
  # Compare alleles
  len1 <- nchar(all1)
  len2 <- nchar(all2)
  if (len1 == 1 && len2 == 1) {
    type = "snp"
  } else if (len1 == len2) {
    type = 'mnp'
  } else if (len1 == 1 || len2 == 1) {
    type = 'indel'
  } else {
    type = 'uneqsub'
  }
  return (list(VarType=type, RefBase=all1, AltBase=all2, RefLength=len1, AltLength=len2))
}

# Redo allele columns in case some have been lost
allele_cols <- colnames(byvar)
allele_cols <- allele_cols[2:length(allele_cols)]


byvar[, {
  res <- class_var(do.call(c, .SD))  # expects list(type=..., all1=..., all2=...)
  for (nm in names(res)) set(byvar, i = .I, j = nm, value = res[[nm]])
  NULL
}, .SDcols = allele_cols, by = .I]

byvar <- byvar[,c("VarID","VarType","RefBase","AltBase","RefLength","AltLength")]

fwrite(byvar, file=paste0(out_pref, ".multi_now_bi.new_allele_info.tsv"))
