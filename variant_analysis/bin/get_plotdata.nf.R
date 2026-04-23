#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(purrr))

# Combine test results from all chromosomes and chunks and for each transfected cell species
# Remove variants with warnings in the model. Leave variants with high skew for filtering out at the next stage
# Compare beta values between cell species
# Correct pvalues
# Add variant details (genotypes and variant types)

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

# Code for getting the relevant data
drop_cols <- c('term', 'conf.level', 'conf.low', 'conf.high', 'statistic', 
               'df.error', 'component', 'original_term', 'has.warning')
get_data <- function(spec, filestem) {
  files <- list.files (pattern = paste0(filestem,"[[:digit:]]+.", spec, ".tsv.gz"))
  data <- rbindlist(lapply(files, fread))
  data <- data[term == 'VarAllele' & !is.na(statistic)]
  
  # remove vars with warnings in model fitting
  warned <- data[(has.warning)]
  warned <- warned[,c("VarID", "estimate", "std.error", "p.value")]
  data <- data[(!has.warning)]
  data <- data[, -drop_cols, with = FALSE]
  
  # Rename columns (except VarID)
  old_names <- setdiff(names(data), "VarID")
  new_names <- paste0(spec, "_", old_names)
  setnames(data, old_names, new_names)
  setnames(warned, old_names, new_names, skip_absent = TRUE)
  
  return(list(warned, data))
}

# Load the normalised stats and model results, then get the model results
allChrom<-list()
allWarns <- list()
for (i in seq_along(chr)) {
  chrid <- chr[i]
  chr_scores <- fread(paste(chrid, out_pref, "meanScores.tsv.gz", sep="."), na.strings = "")
  chr_filestem <- paste0('modelResults.', chrid, "_chunk")
  
  dataList <- lapply(specs, get_data, filestem=chr_filestem)
  warnList <- lapply(dataList, `[[`, 1)
  ChromWarnings <- reduce(warnList, merge, by="VarID", all=TRUE)
  allWarns[[i]] <- ChromWarnings
  
  oneChromList <- lapply(dataList, `[[`, 2)
  rm(dataList)
  oneChrom <- reduce(oneChromList, merge, by="VarID", all=TRUE)
  rm(oneChromList)
  
  # add scores and variant info
  oneChrom <- merge(oneChrom, chr_scores, by = "VarID", all.x = T)
  chr_varfile <- paste0(chrid, ".", out_pref, ".var_details.tsv")
  chr_vardata <- fread(chr_varfile)
  chr_vardata <- chr_vardata[,VarID:AltLengths]
  oneChrom <- merge(oneChrom, chr_vardata, by = c('VarID', 'Chr'), all.x = T)
  chr_typefile <- paste0(chrid, ".", out_pref, ".vartypes.tsv")
  chr_typedata <- fread(chr_typefile)
  oneChrom <- merge(oneChrom, chr_typedata, by = 'VarID', all.x = T)
  
  allChrom[[i]] <- oneChrom
  rm(oneChrom)
}

all_vars<-rbindlist(allChrom)
rm(allChrom)
all_warning <- rbindlist(allWarns)
fwrite(all_warning, file=paste0(out_pref, '.model_warnings.tsv.gz'), sep="\t")

# Do p-value correction for each column
for (species in specs) {
  pval_orig <- paste0(species, "_p.value")
  pval_correct <- paste0(species, "_p.adjust")
  all_vars[, (pval_correct) := p.adjust(get(pval_orig), method = "BH")]
}

# If an interaction term, calculate this using James' beta comparison approach
if (!is.na(interacts[1])) {
  est1 <- paste0(interacts[1], "_estimate")
  est2 <- paste0(interacts[2], "_estimate")
  se1 <- paste0(interacts[1], "_std.error")
  se2 <- paste0(interacts[2], "_std.error")
  all_vars[,interact_z := (
    (get(est1) - get(est2)) / 
      sqrt( (get(se1)*get(se1)) + (get(se2)*get(se2)) )
  )]
  all_vars[, interact_p.value := 2*pnorm(-abs(interact_z))]
  all_vars[, interact_p.adjust := p.adjust(interact_p.value, method="BH")]
}

# Add base position and alleles and update originally-multi variants
# Ref base/length and alt base/length are given but will not actually be ref/alt but instead the lower and higher numbered alt alleles
infocols <- c("VarID", "VarType", "RefLength", "AltLengths", libs)
more_vardata <- fread(paste0(out_pref, ".", "alleles.tsv"), col.names = c("VarID", "Pos","RefBase","AltBase"))
more_vardata <- merge(more_vardata, all_vars[, ..infocols], by = "VarID")
multi_orig <- more_vardata[VarType == 'multi']
multi_recoded <- fread(paste0(out_pref, ".multi_now_bi.new_allele_info.tsv"))
multi <- merge(multi_orig[, .(VarID)], multi_recoded, by = "VarID", all.x = T)
multi_gt <- fread(paste0(out_pref, ".multi_now_bi.recoded_gt.tsv"))
multi <- merge(multi, multi_gt, by = "VarID", all.x = T)
multi[, Pos := NA_integer_]

more_vardata <- more_vardata[VarType != 'multi']
more_vardata[, AltLength := as.numeric(AltLengths)]
more_vardata[, AltLengths := NULL]
more_vardata <- rbind(more_vardata, multi)

infocols2 <- c("VarType", "RefLength", "AltLengths", libs)
all_vars[, (infocols2) := NULL]
all_vars <- merge(all_vars, more_vardata, by="VarID", all.x = T)
all_vars <- setcolorder(all_vars, neworder = c("VarID", "Chr","Pos", "VarType"))

check <- all_vars[is.na(VarType)]
if (nrow(check)) {
  message <- "Some variants were tested but don't have the expected details"
  stop(message)
}


# Split off SNPs and output with allele bases
snps <- all_vars[VarType == 'snp']
snps <- snps[, -c("VarType", "RefLength", "AltLength"), with=FALSE]
fwrite(snps, file=paste0(out_pref, ".snps.tsv.gz"), sep="\t")

# Output all other variants (no allele sequences)
others <- all_vars[VarType != 'snp']
fwrite(others, file=paste0(out_pref, ".non_snps.tsv.gz"), sep="\t", na = "NA")


