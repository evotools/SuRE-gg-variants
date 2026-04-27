#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)  # genomic arithmetic
library(rtracklayer)
library(data.table)
library(BiocParallel)

# Script to combine different variant types and annotate variants with the nearest TSS and PROCap peak
# Do for all tested variants and output as one file labelled with var type
# If ref and alt alleles for indels are >50 bp set to 'Long'
# For non-biallelic this is the total length of the bases string not each allele
# For non-biallelics some columns will be set to NA for compatability.
# Some of this information is available elsewhere (allele lengths, scores, counts)
# Some is not (model estimate and stderr, interaction test)
# Drop information on human+cow combined data

nf_results_path="/exports/eddie/scratch/lplender/SuRE/pipeline_copy2/SuRE_SV/analysis/analysis_nextflow/results/final_stats/"
res1 <- paste0(nf_results_path, "all_cattle.snps.tsv.gz")
res2 <- paste0(nf_results_path, "all_cattle.non_snps.tsv.gz")
multi_file <- paste0(nf_results_path, "all_cattle.multiallelic.tsv.gz")
outdir <- "/home/lplender/scratch/SuRE/analysis_results/sharing_data"

cols_to_drop <- c("Cow_std.error", "Hum_std.error", "interact_z", "all_estimate", "all_std.error","all_p.value","all_p.adjust", "all_0", "all_1", "all_skew")

snps <- fread(res1)
snps <- snps[,VarType := 'snp']
snps <- snps[,RefLength := 1]
snps <- snps[,AltLength := 1]
snps <- snps[, !..cols_to_drop]

nonsnps <- fread(res2)
nonsnps <- nonsnps[, !..cols_to_drop]
nonsnps <- nonsnps[order(Chr, Pos)]
nonsnps <- nonsnps[,RefBase := if_else(nchar(RefBase) <=50, RefBase, "Long")]
nonsnps <- nonsnps[,AltBase := if_else(nchar(AltBase) <=50, AltBase, "Long")]

dat <- rbindlist(list(snps, nonsnps), use.names = T)
rm(snps)
rm(nonsnps)

# fix counts columns (all 4x correct value) - REMOVE IF UPSTREAM CORRECTION!
dat <- dat[,count_0 := (count_0)/4]
dat <- dat[,count_1 := (count_1)/4]

multi <- fread(multi_file)
multi_dropcols <- setdiff(colnames(multi), colnames(dat))
multi_dropcols <- c(multi_dropcols, "count_0","count_1","Cow_0","Cow_1","Hum_0","Hum_1")
multi <- multi[, !..multi_dropcols]
multi <- multi[,AltBase := if_else(nchar(AltBase) <=50, AltBase, "Long")] #combined not individual lengths (plus comma)
multi <- multi[,VarType := 'multi']

dat <- rbindlist(list(dat, multi), use.names = T, fill = T)
dat <- as_tibble(dat)

fdrThresh <- 0.05
skew_thresh <- 20

dat <- dat %>%
  rename_with(tolower) %>%
  rename(refbase = "ref", altbase = "alt") %>%
  # Add significance flags for each test
  mutate(
    cow_signif = if_else(cow_p.adjust <= fdrThresh, "Y", "N"),
    hum_signif = if_else(hum_p.adjust <= fdrThresh, "Y", "N"),
    interact_signif = if_else(interact_p.adjust <= fdrThresh, "Y", "N")
  ) %>% 
  filter(!is.na(pos))

# Remove variants with high skew (to avoid tripping up later)
# and those that were not tested in cow
# Keep high human skew for now, remove at the end for second output
dat <- dat %>% filter(cow_skew<=skew_thresh, !is.na(cow_p.adjust))

message ("Writing significant variants outputs")
cow_signif <- dat %>% filter(cow_signif == "Y", vartype != 'multi') %>% 
  select(varid, chr, pos, vartype, ref, alt, cow_p.adjust, hum_p.adjust, hum_skew) %>%
  arrange(chr, pos)
fwrite(cow_signif, file=paste0(outdir, "cow_signif.not_hum_skew_flt.tsv.gz"), sep="\t", na = 'NA', quote = F)

hum_signif <- dat %>% filter(hum_signif == "Y", vartype != 'multi', hum_skew <= skew_thresh) %>% 
  select(varid, chr, pos, vartype, ref, alt, hum_p.adjust, cow_p.adjust, hum_skew) %>%
  arrange(chr, pos)
fwrite(hum_signif, file=paste0(outdir, "hum_signif.tsv.gz"), sep="\t", na = 'NA', quote = F)

xs_signif <- dat %>% filter(interact_signif == "Y", vartype != 'multi', hum_skew <= skew_thresh) %>%
  select(varid, chr, pos, ref, alt, cow_p.adjust, hum_p.adjust, interact_p.adjust)
fwrite(xs_signif, file=paste0(outdir, "xs_signif.tsv.gz"), sep="\t", na = 'NA', quote = F)


message("Adding TSS information")
## -----------------------------------------------------------
# Add TSS information to the data
## -----------------------------------------------------------
gtf_url  <- "https://ftp.ensembl.org/pub/release-115/gtf/bos_taurus/Bos_taurus.ARS-UCD2.0.115.gtf.gz"
## import() understands remote .gz files – no explicit download needed
gtf <- rtracklayer::import(gtf_url)
## keep one record per *transcript* (or per *gene* – your choice)
tx   <- gtf[gtf$type == "transcript"]
## TSS = first base of the transcript respecting strand
tss  <- resize(tx, width = 1, fix = "start")   # start on "+" , end on "-"
## we only need chr + 1-bp position + gene info
tss  <- unique(tss)           # many transcripts share the same TSS

# edit here to allow longer than 1 bp vars
gr_var <- GRanges(seqnames = dat$chr,
                  ranges   = IRanges(dat$pos, dat$pos+dat$reflength-1),
                  strand   = "*")          # variants are un-stranded

hits <- distanceToNearest(gr_var, tss, ignore.strand = TRUE)

## absolute distance returned by distanceToNearest()
abs_d <- mcols(hits)$distance

## position of variant and the matched TSS
var_pos <- start(gr_var)[queryHits(hits)]
tss_pos <- start(tss)[subjectHits(hits)]

## signed distance:  +ve if the TSS is downstream of the variant
signed_d <- tss_pos - var_pos      # downstream → positive, upstream → negative                  

## Extract gene information from the nearest TSS
## Note: subjectHits(hits) gives the index of the nearest TSS in the tss GRanges
nearest_tss_indices <- subjectHits(hits)
gene_id <- mcols(tss)$gene_id[nearest_tss_indices]
gene_biotype <- mcols(tss)$gene_biotype[nearest_tss_indices]
gene_name <- mcols(tss)$gene_name[nearest_tss_indices]

dat <- dat %>%
  mutate(nearest_tss_chr = as.character(seqnames(tss)[nearest_tss_indices]),
         nearest_tss_pos = tss_pos,
         dist_to_tss     = signed_d,
         nearest_gene_id = gene_id,
         nearest_gene_biotype = gene_biotype,
         nearest_gene_name = gene_name)

## -----------------------------------------------------------
# Add PROCap information to the data
## -----------------------------------------------------------

message("Adding PROcap information")
# Windows of 500 bp each side of the variant -- flanks for >1 bp
win_gr <- copy(gr_var)
start(win_gr) <- start(win_gr) - 500
end(win_gr) <- end(win_gr) + 500
win_gr <- trim(win_gr)

# -------------------------------------------------------------------
# Extract max PROCap score in a 1-kb window centred on each variant.
# This is the slowest step of the script (~20 min with 8 cores).
#   To accelerate, we configure a parallel backend using BiocParallel.
#   On Windows, SnowParam is used; on Linux/macOS, use MulticoreParam.
#   Not sure this is doing much but doesnt hurt.
# -------------------------------------------------------------------
if (.Platform$OS.type == "unix") {
  bpparam <- MulticoreParam(workers = parallel::detectCores() - 1)
} else {
  bpparam <- SnowParam(workers = parallel::detectCores() - 1, type = "SOCK")
}

bw_file <- "/home/lplender/scratch/SuRE/analysis_results/data/merged_RPI2_BOAEC-D1_RPI5.deduplicated_without_bl.bw"

bw <- BigWigFile(bw_file)
## get the maximum score in every window
message("starting PROcap lookup")
start_time <- Sys.time()
summ <- summary(bw, win_gr, type = "max", bpparam = bpparam)   # returns DataFrame with $score
message(sprintf("Run time: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
scores <- unlist(summ)$score

dat$PROCap_max <- unlist(summ)$score        # NA means no data → treat as 0 if you like


message("Adding frequency information")
## Add pangenie allele frequencies
pg_freqsfiles <- list.files(path="/exports/eddie/scratch/lplender/SuRE/pangenie_frequencies/", full.names = T)

names(pg_freqsfiles) <- c("Taurine","Indicine")
col_nms_pg <- c("chr","pos","varid","AF","AN")
pg_freqs <- data.table::rbindlist(lapply(pg_freqsfiles, data.table::fread, col.names = col_nms_pg, na.strings = "."), idcol = "breed")
pg_noid <- pg_freqs[is.na(varid)]
pg_freqs <- pg_freqs[!is.na(varid)]
pg_freqs <- dcast(pg_freqs, formula = varid + chr + pos ~ breed, value.var = c("AF","AN"))
pg_freqs[, chr := as.numeric(chr)]
pg_freqs <- pg_freqs[!is.na(chr)]

dat <- left_join(dat, pg_freqs)

message("Writing main outputs")

# Write output
outfile1 <- paste0(outdir, "/all_vars.annotated.cow_skew_filtered.tsv.gz")
fwrite(dat, outfile1, quote = F, sep="\t", na="NA")

# Remove variants with high skew in human
dat <- dat %>% filter(hum_skew<=skew_thresh, !is.na(hum_p.value))
outfile2 <- paste0(outdir, "/all_vars.annotated.both_skew_filtered.tsv.gz")
fwrite(dat, outfile2, quote = F, sep="\t", na="NA")

