library(data.table)
library(stringr)

# MPRA results come from Annogen in a large table with only a few columns we need for this analysis
# Extract the fragment and variant details, input count, cDNA counts
# For this experiment, 'input' is the pDNA count not the iPCR count. This is summed across the two replicates

# Assumes only one alt allele per variant, will need to update if this changes

args <- commandArgs(trailingOnly = TRUE)
config <- args[1]
source(config)

outfile <- ".Mamm.tsv.gz"
keepcols <- c("Element", "Strand", "rsID", "refalt", "Barcode", pdna_cols)

# Get relevant columns from each Barcodes.txt.gz file
biorep1 <- fread(infiles_batch1[1])
keep1 <- c(keepcols, intersect(colnames(biorep1), reps))
biorep1 <- biorep1[,..keep1]
biorep2 <- fread(infiles_batch1[2])
keep2 <- c(keepcols, intersect(colnames(biorep2), reps))
biorep2 <- biorep2[,..keep2]

# Get the PAEP fragments from the second batch
paeprep1 <- fread(infiles_batch2[1])
paeprep1 <- paeprep1[str_detect(Element, 'paep')][,..keep1]
paeprep2 <- fread(infiles_batch2[2])
paeprep2 <- paeprep2[str_detect(Element, 'paep')][,..keep2]

biorep1 <- rbind(biorep1, paeprep1)
biorep2 <- rbind(biorep2, paeprep2)

# Merge the bio replicates and check non-cDNA columns agree
if (nrow(biorep1) != nrow(biorep2)) {
  message("Unequal replicates")
}
countstable <- merge(biorep1, biorep2, by=keepcols)
if (nrow(biorep1) != nrow(countstable)) {
  message("Replicates don't match")
}

# This script is only for biallelic vars, check this is the case here (only one 'alt' string)
alts <- countstable[str_detect(string=refalt, pattern='alt')]$refalt %>% unique()
if (!length(alts == 1)) {
  message("More than one alt allele! Don't use")
}

# Sum then drop pDNA columns (if >1)
countstable[, input := rowSums(.SD), .SDcols = pdna_cols]
countstable[, (pdna_cols) := NULL]

# Sum the cDNA columns for each bio rep
countstable[, Rep1 := rowSums(.SD), .SDcols = reps[1:2]]
countstable[, Rep2 := rowSums(.SD), .SDcols = reps[3:4]]
countstable[, (reps) := NULL]


# Add VarAllele column with Ref/Alt converted to 0/1 to follow previous model format
# Since each alt will be tested separately with the reference, all alts can be converted to 0
countstable <- countstable[, VarAllele := fifelse(refalt=="ref", "0", "1")]

# Add extra columns needed by the next script
countstable[, cellSpecies := 'Cow'][, library := 'n']

# Output with strands combined (use rsID as-is)
strand_comb <- copy(countstable)[, c("Element", "Strand","refalt") := NULL]
setnames(strand_comb, old="rsID", new="VarID")
comb_dat <- melt(
  strand_comb,
  id.vars = c("VarID","VarAllele", "input", "library","Barcode","cellSpecies"), # keep these
  measure.vars = c("Rep1","Rep2"), # columns to pivot
  variable.name = "Rep",
  value.name = "cDNA", # new column holding the numeric count
  variable.factor = FALSE # keep 'cells' as plain character
)

fwrite(comb_dat, paste0(data_outdir, "nonstranded", outfile), sep="\t")

# Stranded output
stranded <- copy(countstable)[, c("Element","refalt") := NULL]
stranded[, VarID := paste(rsID, Strand, sep="_")][, c("rsID","Strand") := NULL]
stranded_dat <- melt(
  stranded,
  id.vars = c("VarID","VarAllele", "input", "library","Barcode","cellSpecies"), # keep these
  measure.vars = c("Rep1","Rep2"), # columns to pivot
  variable.name = "Rep",
  value.name = "cDNA", # new column holding the numeric count
  variable.factor = FALSE # keep 'cells' as plain character
)

fwrite(stranded_dat, paste0(data_outdir, "stranded", outfile), sep="\t")
