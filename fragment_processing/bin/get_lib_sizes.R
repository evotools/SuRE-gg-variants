#!/usr/bin/env Rscript

library(data.table)

# Process counts table files to get lib size files
# Start from the iPCR column and do all to the end

args <- commandArgs(trailingOnly = TRUE)
input.file <- args[1]
output.file <- args[2]

# Read the data
data <- fread(input.file)

# Work out which columns to get
columns <- colnames(data)
ipcr_col <- which(columns=='iPCR')
to_sum <- columns[ipcr_col:length(columns)]

# Calculate column totals for numeric columns only
column_totals <- data[, lapply(.SD, function(x) if (is.numeric(x)) sum(x, na.rm=TRUE) else NA), .SDcols = to_sum]

# Add a column with the number of rows in data
column_totals[, unique_frags := nrow(data)]

# Write output
fwrite(column_totals, file=output.file, sep="\t", quote=FALSE)
