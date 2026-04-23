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
chrid <- as.character(args[3])
chunkfile <- args[4]
out_pref <- args[5]
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
j <- sub(".*_(\\d+)\\.tsv.gz$", "\\1", chunkfile)

by_var_flat <- fread(file=chunkfile)
by_var_flat[, ':=' (library = NULL, Barcode = NULL) ]

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
rm(by_var_flat)
gc()
var_data <- tibble::tibble(VarID = names(var_data), data  = purrr::map(unname(var_data), tibble::as_tibble))

# Cleaned up code to set model and get data
model_fml <- if (all) {
  cDNA ~ offset(log(input)) + cellSpecies + VarAllele + cellSpecies:VarAllele | Rep
} else {
  cDNA ~ offset(log(input)) + Rep + VarAllele | 1
}

safe_quiet_model <- safely(quietly(function(df) {
  tidy_zeroinfl(
    zeroinfl(model_fml, data = df, dist = "negbin")
  )
}))


# Extract list of data frames
data_list <- var_data$data

message(sprintf("Got data, starting to run models: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

# Parallel processing
t_model <- system.time({safe_models <- mclapply(data_list, safe_quiet_model, mc.cores = n_cores)})
cat(sprintf("Model processing time: %7.2f seconds  (workers = %d)\n", t_model["elapsed"], n_cores))

message(sprintf("Finished models: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
# Add results back to var_data
var_data$safe.model <- safe_models

# Finish the processing including adding annotation for failed/warnings
tests <- var_data %>%
  mutate(
    failed = map_lgl(safe.model, ~ !is.null(.x$error)),
    model = map(safe.model, ~ .x$result$result),
    has.warning = map_lgl(safe.model, ~ {
      w <- .x$result$warnings
      length(w) > 0
    })
  ) %>%
  select(-safe.model, -data)

tests.failed <- tests %>% filter(failed) %>% pull(VarID)
tests.failed <- if (length(tests.failed) == 0) NA else tests.failed
write(tests.failed, file=paste0("model_failed.", chrid, "_chunk", j, ".", spec, ".tsv"))

cat("Failed VarIDs:\n")
cat(tests.failed, sep = "\n")

tests.ok <- tests %>% filter(!failed) %>% select(-failed)
mod_res<- tests.ok %>% unnest(model)

write_tsv(mod_res, paste0("modelResults.", chrid, "_chunk", j, ".", spec, ".tsv.gz"))


## ── 9.  Clean-up ──────────────────────────────────────────────
future::plan(sequential)
message("Finished chunk chr", chrid, "-", j, " (", spec, ")")