
library(MASS)   # for glm.nb()
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(pscl))
library(nplyr)
library(parallel)
library(furrr)
library(future)
library(broom)
library(broom.helpers)

# Run the allele-specific expression model on MPRA fragment data, labelled by strand
# Models:
# 1) All fragments, modelled with Strand as covariate (where data for both strand)
# 2) Both strands separately

args <- commandArgs(trailingOnly = TRUE)
chunkfile <- args[1]
outprefix <- args[2]
config <- args[3]
source(config)

## ── 1.  Cores & parallel plan (forked) ────────────────────────
# 1a. how many cores did the scheduler give us?
n_cores <- as.integer(Sys.getenv("NSLOTS", "1"))
if (is.na(n_cores) || n_cores < 1)
  n_cores <- future::availableCores()
# 1b. do not exceed the number of physical cores
n_cores <- min(n_cores, parallel::detectCores(logical = FALSE))
message("Using ", n_cores, " forked R workers")
# 1c. multicore (fork) if supported, otherwise fall back gracefully
if (future::supportsMulticore()) {
  future::plan(multicore, workers = n_cores)
} else {
  warning("Multicore futures not supported on this platform; using sequential")
  future::plan(sequential)
}

data.table::setDTthreads(n_cores)        # keep data.table aligned
options(future.rng.onMisuse = "ignore")  # suppress RNG warning

# Set models
nb_model_fml <- cDNA ~ offset(log(input)) + Rep + VarAllele
safe_quiet_model_nb <- safely(quietly(function(df) {
  tidy(
    glm.nb(nb_model_fml, data = df)
  )
}))

nb_model_fml_strand <- cDNA ~ offset(log(input)) + Rep + VarAllele + Strand
safe_quiet_model_nbs <- safely(quietly(function(df) {
  tidy(
    glm.nb(nb_model_fml_strand, data = df)
  )
}))


# Actual data handling
var_data_flat <- fread(file=chunkfile)
var_data_flat[, ':=' (library = NULL, Barcode = NULL) ]

# Check fragments for both alleles
allele_counts <- as_tibble(var_data_flat) %>% group_by(VarID) %>%
  count(VarAllele) %>% count(VarID)
only1 <- allele_counts %>% filter(n<2)
if (nrow(only1) > 0) {
  message ("Some variants have only one allele, removing")
  write(only1$VarID, file=paste0(outprefix, "oneAll.txt"))
  var_data_flat <- var_data_flat[!(VarID %in% only1$VarID)]
}

# Get data for strand-separated models
var_data_stranded <- split(var_data_flat, by = "VarID", keep.by = FALSE, flatten = FALSE)
var_data_stranded <- tibble::tibble(VarID = names(var_data_stranded), data  = purrr::map(unname(var_data_stranded), tibble::as_tibble))

# add strand as column
var_data_flat <- var_data_flat[, c("VarID", "Strand") := tstrsplit(VarID, "_", fixed = TRUE, keep = c(1L, 2L))]
var_data <- split(var_data_flat, by = "VarID", keep.by = FALSE, flatten = FALSE)
var_data <- tibble::tibble(VarID = names(var_data), data  = purrr::map(unname(var_data), tibble::as_tibble))

# Count strand. If only 1 then strand model will not run properly so annotate
vars_strandcount <- as_tibble(var_data_flat) %>% distinct(VarID, Strand) %>% 
  count(VarID) %>%
  mutate(withStrand = ifelse(n == 1, "N", ifelse(n==2, "Y", "??")), .keep='unused')
var_data <- left_join(var_data, vars_strandcount)



one_model <- function(df, thismodel, failprefix) {
  # Extract list of data frames
  data_list <- df$data
  # Run the model
  t_model <- system.time({safe_models <- mclapply(data_list, thismodel, mc.cores = n_cores)})
  message(cat(sprintf("Model processing time: %7.2f seconds  (workers = %d)\n", t_model["elapsed"], n_cores)))
  # Add back to var data
  df$safe.model <- safe_models
  
  # Finish the processing including adding annotation for failed/warnings
  tests <- df %>%
    mutate(
      failed = map_lgl(safe.model, ~ !is.null(.x$error)),
      model = map(safe.model, ~ .x$result$result),
      has.warning = map_lgl(safe.model, ~ {
        w <- .x$result$warnings
        length(w) > 0
      })
    ) %>% 
    dplyr::select(-safe.model, -data)
  tests.failed <- tests %>% filter(failed) %>% pull(VarID)
  tests.failed <- if (length(tests.failed) == 0) NA else tests.failed.nbs
  write(tests.failed, file=paste0(outprefix, failprefix, "model_failed.tsv"))
  
  # Get the results for VarAllele
  tests.ok <- tests %>% filter(!failed) %>% dplyr::select(-failed)
  mod_res <- tests.ok %>% unnest(model) %>% filter(term=='VarAllele')
  return(mod_res)
}

# For combined strand data, run model with strand covariate if both strands present, otherwise just for the one that's there
# P-value correction (process data_subset, below) will use the p-value from the appropriate test for each variant
bothstrands <- var_data %>% filter(withStrand == "Y")
bothstrands_res <- one_model(bothstrands, safe_quiet_model_nbs, "both_strands.")
onestrand <- var_data %>% filter(withStrand == "N")
onestrand_res <- one_model(onestrand, safe_quiet_model_nb, "one_strand.")

# Merge the model results
keepcols <- c("VarID", "estimate", "p.value", "has.warning", "withStrand")
comb_results <- bind_rows(bothstrands_res %>% dplyr::select(all_of(keepcols)), onestrand_res %>% dplyr::select(all_of(keepcols)))

# Get the relevant rsIDs
elements_info <- read_tsv(elements_infofiles, show_col_types = FALSE) %>% distinct(rsID, Element, Strand)
comb_results <- left_join(comb_results, elements_info %>% distinct(rsID, Element), by=c("VarID"="rsID"))
comb_results <- comb_results %>% relocate(Element, .after = VarID)

# Do a separate model for each strand
sep_strand_res <- one_model(var_data_stranded, safe_quiet_model_nb, "sep_strands.")
sep_strand_res <- sep_strand_res %>% separate_wider_delim(cols=VarID, delim = "_", names = c("VarID","Strand"))
keepcols <- c("VarID", "Strand","estimate", "p.value", "has.warning")
sep_strand_res <- left_join(sep_strand_res %>% dplyr::select(all_of(keepcols)), elements_info %>% distinct(rsID, Element), by=c("VarID"="rsID")) %>%
  relocate(Element, .after = VarID)

# Do p-value correction for each set of variants then output
data_subset <- function(df, string) {
  sub_res <- df %>% filter(str_detect(Element, string))
  sub_res <- sub_res %>% mutate(
    p.adjust = p.adjust(p.value, method="BH"),
    signif = ifelse(p.adjust <= 0.05, "Y","N")
  )
  return(sub_res)
}

comb_list <- list()
sep_list <- list()
for (i in seq_along(subset_tags)) {
  results <- data_subset(comb_results, subset_tags[i])
  write_tsv(results, paste0(outprefix, "modelRes.", subset_names[i], ".tsv"))
  results <- results %>% mutate(subset=subset_tags[i])
  comb_list[[i]] <- results
  
  sep_strand_results <- data_subset(sep_strand_res, subset_tags[i])
  write_tsv(results, paste0(outprefix, "modelRes.sep_strands.", subset_names[i], ".tsv"))
  sep_strand_results <- sep_strand_results %>% mutate(subset=subset_tags[i])
  sep_list[[i]] <- sep_strand_results
}

comb_tagged <- rbindlist(comb_list)
write_tsv(comb_tagged, paste0(outprefix, "strand_combined.all_subsets.tsv"))
comb_counts <- comb_tagged %>% count(signif, subset) %>%
  pivot_wider(id_cols=subset, names_from = signif, values_from = n) %>% 
  mutate(total=Y+N)
write_tsv(comb_counts, paste0(outprefix, "comb_counts.tsv"))

sep_tagged <- rbindlist(sep_list)
write_tsv(sep_tagged, paste0(outprefix, "strand_sep.all_subsets.tsv"))
sep_counts <- sep_tagged %>% count(signif, subset, Strand)  %>%
  pivot_wider(id_cols=c(subset, Strand), names_from = signif, values_from = n) %>% 
  mutate(total=Y+N)
write_tsv(sep_counts, paste0(outprefix, "sep_strand_counts.tsv"))

# Counts for whether variant is signif on both strands, one or neither
sep_tagged_wide <- sep_tagged %>% 
  pivot_wider(id_cols = c(VarID, Element, subset), names_from=Strand, values_from=signif)
sep_wide_counts <- sep_tagged_wide %>% count(subset,`+`,`-`) %>%
  replace_na(list(`-`="0",`+`="0")) %>%
  mutate(class = case_when(
    `+` == "Y" & `-` == "N" ~ '+only',
    `+` == "N" & `-` == "Y" ~ '-only',
    `+` == "Y" & `-` == "Y" ~ 'both',
    `+` == "N" & `-` == "N" ~ 'neither',
    `+` == "0" & `-` == "N" ~ 'neither,no+',
    `+` == "0" & `-` == "Y" ~ '-only,no+',
    `+` == "N" & `-` == "0" ~ 'neither,no-',
    `+` == "Y" & `-` == "0" ~ '+only,no-'
  )) %>%
  pivot_wider(id_cols=subset, names_from = class, values_from=n)
write_tsv(sep_wide_counts, paste0(outprefix, "sep_strand_counts.by_var.tsv"))
