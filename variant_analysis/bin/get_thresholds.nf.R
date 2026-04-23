#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))

# Calculate thresholds for variant filtering using data combined from all chromosomes
# Also plot distributions and thresholds for checking
# Thresholds are for total iPCR coverage, proportion of reference allele (non-imputed only) and proportion of imputed (fragments with any imputated only), calculated for each library_genotype

# Stop on any error immediately
options(error = function() {
  cat("Error occurred at:", Sys.time(), "\n")
  cat("Error message:", geterrmessage(), "\n")
  traceback()
  quit(status = 1)  # Exit with error code
})
#options(warn = 2)  # Don't turn warnings into errors!

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
out_pref <- args[2]
source(config_file)

files <- list.files(pattern="counts.rds")

# Read and combine files from all chr
load_rdata <- function(file) {
  fragRatios <- readRDS(file)
  chrid <- str_extract(file, "^[^.]+")
  fragRatios <- fragRatios %>% mutate(chr = chrid)
  return(fragRatios)
}
allcounts <- list_rbind(map(files, load_rdata)) %>% ungroup
#allcounts.cp <- allcounts

# Calculate thresholds for each library+genotype combination. Calculate genome-wide for filtering and per-chr for checking

pdf(file=paste0(out_pref,".threshold_plots.pdf"), width=11)

# For data sets with some imputed variants
if ("impProp" %in% colnames(allcounts)) {
  # First of all calculate thresholds for all imputed proportions (across the whole variant). Only want upper threshold as fine if no fragments are imputed.
  threshImp <- allcounts %>% group_by(genotype, library) %>% 
    summarise(meanImp=mean(impProp, na.rm=TRUE), upImp=meanImp+2.58*sd(impProp, na.rm=TRUE))
  threshImp.bychr <- allcounts %>% group_by(genotype, library, chr) %>% 
    summarise(meanImp=mean(impProp, na.rm=TRUE), upImp=meanImp+2.58*sd(impProp, na.rm=TRUE), n=n())
  plot1 <- ggplot(allcounts, mapping=aes(x=impProp)) + 
    geom_density() + 
    geom_vline(data = threshImp, aes(xintercept = meanImp), linetype='dashed') +
    geom_vline(data = threshImp, aes(xintercept = upImp), linetype='dashed', colour='red') +
    facet_wrap(vars(library, genotype), scales='free') +
    labs(title = "Proportion of fragments that are imputed")
  print(plot1)
  
  # Output
  write_tsv(threshImp, file=paste0(out_pref,".thresholds.imp.tsv"))
  write_tsv(threshImp.bychr, file=paste0(out_pref,".by_chr_thresholds.imp.tsv"))
}

# Mean and mean+/- SD for reference proportion (fraction of all fragments assigned to a variant+library+genotype that are 0)
threshRatio<-allcounts %>% group_by(genotype, library) %>% 
  summarise(meanRatio=mean(propRef, na.rm=TRUE), downRatio=meanRatio-2.58*sd(propRef, na.rm=TRUE), upRatio=meanRatio+2.58*sd(propRef, na.rm=TRUE))
threshRatio.bychr <- allcounts %>% group_by(genotype, library, chr) %>% 
  summarise(meanRatio=mean(propRef, na.rm=TRUE), downRatio=meanRatio-2.58*sd(propRef, na.rm=TRUE), upRatio=meanRatio+2.58*sd(propRef, na.rm=TRUE), n=n())

plot2 <- ggplot(allcounts, mapping=aes(x=propRef)) + geom_density() + 
  geom_vline(data = threshRatio, aes(xintercept = meanRatio), linetype='dashed') +
  geom_vline(data = threshRatio, aes(xintercept = downRatio), linetype='dashed', colour='blue') +
  geom_vline(data = threshRatio, aes(xintercept = upRatio), linetype='dashed', colour='red') +
  facet_wrap(vars(library, genotype), scales='free') +
  labs(title = "Proportion of non-imputed fragments that are ref (0)")
print(plot2)

# Mean and mean+/- SD for total iPCR count. Hard filter first for >=1000
threshCov <- allcounts %>% filter(varInput < 1000) %>% 
  group_by(genotype, library) %>% 
  summarise(meanCov=mean(varInput), downCov=meanCov-2.58*sd(varInput), upCov = meanCov+2.58*sd(varInput))
threshCov.bychr <- allcounts %>% ungroup %>% filter(varInput < 1000) %>%
  group_by(genotype, library,chr) %>% 
  summarise(meanCov=mean(varInput), downCov=meanCov-2.58*sd(varInput), upCov = meanCov+2.58*sd(varInput))

plot3 <- ggplot(allcounts %>% filter(varInput < 1000), mapping=aes(x=varInput)) + geom_density() + 
  geom_vline(data = threshCov, aes(xintercept = meanCov), linetype='dashed') +
  geom_vline(data = threshCov, aes(xintercept = downCov), linetype='dashed', colour='blue') +
  geom_vline(data = threshCov, aes(xintercept = upCov), linetype='dashed', colour='red') +
  facet_wrap(vars(library, genotype), scales='free') +
  labs(title = "Sum of iPCR counts for all fragments")
print(plot3)

# Output
write_tsv(threshRatio, file=paste0(out_pref,".thresholds.ratio.tsv"))
write_tsv(threshRatio.bychr, file=paste0(out_pref,".by_chr_thresholds.ratio.tsv"))
write_tsv(threshCov, file=paste0(out_pref,".thresholds.iPCR.tsv"))
write_tsv(threshCov.bychr, file=paste0(out_pref,".by_chr_thresholds.iPCR.tsv"))

dev.off()


