#!/usr/bin/env Rscript

start_time <- Sys.time()

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(pscl))
suppressMessages(library(nplyr))
suppressMessages(library(parallel))
suppressMessages(library(furrr))
suppressMessages(library(future))
suppressMessages(library(broom))
suppressMessages(library(broom.helpers))
suppressMessages(library(lmtest))

# Stop on any error immediately
options(error = function() {
  cat("Error occurred at:", Sys.time(), "\n")
  cat("Error message:", geterrmessage(), "\n")
  traceback()
  quit(status = 1)  # Exit with error code
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
spec<- args[2]
chunkfile <- args[3]
out_pref <- args[4]
source(config_file)

## ── 1.  Cores & parallel plan (forked) ────────────────────────
# Copied from James script runModelChunks.R
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

# Actual data handling
# get chunk number from filename
j <- sub("chunk_multi.(\\d+)\\.tsv.gz$", "\\1", chunkfile)

by_var_flat <- fread(file=chunkfile)
by_var_flat[, ':=' (library = NULL, Barcode = NULL) ]

# Encode allele as factor
by_var_flat[ , VarAllele := as.factor(VarAllele)]

# filter for species before grouping
all <- tolower(spec) == "all"
var_data_flat  <- if (all) {
  by_var_flat
} else {
  by_var_flat[cellSpecies == spec]
}
rm(by_var_flat)
gc()

var_data <- split(var_data_flat, by = "VarID", keep.by = FALSE, flatten = FALSE)
rm(var_data_flat)
gc()
var_data <- tibble::tibble(VarID = names(var_data), data  = purrr::map(unname(var_data), tibble::as_tibble))

# Set model
fml_var <- if (all) {
  cDNA ~ offset(log(input)) + cellSpecies + VarAllele | Rep
} else {
  cDNA ~ offset(log(input)) + Rep + VarAllele | 1
}

fml_novar <- if (all) {
  cDNA ~ offset(log(input)) + cellSpecies | Rep
} else {
  cDNA ~ offset(log(input)) + Rep | 1
}

safe_quiet_model_var <- safely(quietly(function(df) {
  zeroinfl(fml_var, data = df, dist = "negbin")
}))

safe_quiet_model_novar <- safely(quietly(function(df) {
  zeroinfl(fml_novar, data = df, dist = "negbin")
}))

# Extract list of data frames
data_list <- var_data$data
var_data <- var_data %>% select(-data)

message(sprintf("Got data, starting to run models: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Parallel processing
t_model1 <- system.time({var_models <- mclapply(data_list, safe_quiet_model_var, mc.cores = n_cores)})
cat(sprintf("With-var model processing time: %7.2f seconds  (workers = %d)\n", t_model1["elapsed"], n_cores))

# Add results back to var_data
var_data$var_result <- var_models
rm(var_models)
gc()

t_model2 <- system.time({novar_models <- mclapply(data_list, safe_quiet_model_novar, mc.cores = n_cores)})
cat(sprintf("With-var model processing time: %7.2f seconds  (workers = %d)\n", t_model2["elapsed"], n_cores))
var_data$novar_result <- novar_models
rm(novar_models)
gc()

message(sprintf("Finished models: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Remove rows where either model failed
var_data <- var_data %>%
  mutate(
  failed_var = map_lgl(var_result, ~ !is.null(.x$error)),
  failed_novar = map_lgl(novar_result, ~ !is.null(.x$error)),
)

tests_failed <- var_data %>% filter(failed_var | failed_novar) %>% select(VarID, failed_var, failed_novar)
var_data <- var_data %>% filter(!failed_var & !failed_novar) %>%
  select(-failed_var, -failed_novar)

# Extract models, run likelihood ratio test and add annotation for warnings
var_data <- var_data %>% mutate(var_model = map(var_result, ~ .x$result$result), 
                             novar_model = map(novar_result, ~ .x$result$result),
                             var_warnings = map_lgl(var_result, ~ {
                               w <- .x$result$warnings
                               length(w) > 0
                             }),
                             novar_warnings = map_lgl(novar_result, ~ {
                               w <- .x$result$warnings
                               length(w) > 0
                             })) %>%
  select(-var_result, -novar_result) %>%
  mutate(comparison=map2(var_model, novar_model, ~ lrtest(.x, .y))) %>%
  select(-var_model, -novar_model)

mod_res <- var_data %>% unnest(comparison) %>% 
  rename("pval.lrtest"="Pr(>Chisq)") %>% filter(!is.na(pval.lrtest)) %>%
  select(VarID:novar_warnings, pval.lrtest)

tests_failed <- if (length(tests_failed) == 0) NA else tests_failed
write_tsv(tests_failed, file=paste0("model_failed.chunk_multi.", j, ".", spec, ".tsv"))

write_tsv(mod_res, paste0("modelResults.chunk_multi.", j, ".", spec, ".tsv.gz"))


## ── 9.  Clean-up ──────────────────────────────────────────────
future::plan(sequential)
message("Finished chunk ", j, " (", spec, ")")