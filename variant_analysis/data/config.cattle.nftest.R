# Variables for analysis
# chromosome names
#chr <- as.character(c(1:29))
chr <- c("28","29")
# Names of SuRE libraries (as in var counts tables)
libs <- c('Bi-83-A','BOAEC-D1')

# Names of SuRE transfection replicates (should match var counts tables)
reps <- c("BOAEC_D1","BOAEC_D2","HAEC_D1_2","HAEC_D1_3")

# Aliases for SuRE transfection replicates to use in this analysis, in same order as in reps
# Should be in the format species_RepID, with RepID unique for each not just within the species, e.g. Cow_Rep1, Cow_Rep2, Hum_Rep3...
# species should match any abbreviations in the nextflow config file
reps_rename <- c("Cow_Rep1","Cow_Rep2","Hum_Rep3", "Hum_Rep4")

# sets of data to run -- 'all' to do complete data set
# And two of these to test for an interaction between. If doing >1 pair might have to rewrite. If none set as NA
specs <- c('Cow', 'Hum')
interacts <- c('Cow', 'Hum')

# Filtering options
minMAPQ <- 0  # only tested if non-0
filter.cav <- FALSE
cav.excl <- c("I","II","E") # E should already have been excluded anyway

chunk_count <- 5
pval_cutoff <- 0.05

skew_threshold <- 20  # this is calculated for each cell type individually in normalise_SuREdata.R
walks_ids <- c("Bi_83_A", "BOAEC_D1")
