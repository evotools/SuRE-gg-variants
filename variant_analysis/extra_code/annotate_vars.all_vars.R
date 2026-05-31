#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))  # genomic arithmetic
suppressMessages(library(rtracklayer))
suppressMessages(library(data.table))
suppressMessages(library(BiocParallel))

# Script to combine different variant types and annotate variants with the nearest TSS and PROCap peak
# Do for all tested variants and output as one file labelled with var type
# If ref and alt alleles for indels are >50 bp set to '>50'
# For multiallelics, values for all alt alleles 1-4 (whether present or not) are combined into one column for each category (fragment counts, alleles, lengths etc).

# Data required in data_dir:
# PRO-Cap bigwig file
# Alt allele frequencies from pangenie
# Ancestral alleles

args <- commandArgs(trailingOnly = TRUE)
nf_results_path <- args[1]
outdir <- args[2]
data_dir <- args[3]

res1 <- paste0(nf_results_path, "all_cattle.snps.tsv.gz")
res2 <- paste0(nf_results_path, "all_cattle.non_snps.tsv.gz")
multi_file <- paste0(nf_results_path, "all_cattle.multiallelic.tsv.gz")

cols_to_drop <- c("Cow_std.error", "Hum_std.error", "interact_z", "all_estimate", "all_std.error","all_p.value","all_p.adjust", "all_0", "all_1", "all_skew")

message("Loading variant results")
snps <- fread(res1)
snps <- snps[,VarType := 'snp']
snps <- snps[,RefLength := 1]
snps <- snps[,AltLength := 1]
snps <- snps[, !..cols_to_drop]

nonsnps <- fread(res2)
nonsnps <- nonsnps[, !..cols_to_drop]
nonsnps <- nonsnps[order(Chr, Pos)]
nonsnps <- nonsnps[,RefBase := if_else(nchar(RefBase) <=50, RefBase, ">50")]
nonsnps <- nonsnps[,AltBase := if_else(nchar(AltBase) <=50, AltBase, ">50")]

dat_all <- rbindlist(list(snps, nonsnps), use.names = T)

# fix counts columns (all 4x correct value -- in the file I counted rows from there was one row per transfection per fragment) - REMOVE IF UPSTREAM CORRECTION!
dat_all <- dat_all[,count_0 := (count_0)/4]
dat_all <- dat_all[,count_1 := (count_1)/4]

setnames(dat_all, old=c("count_0","count_1", "Cow_0", "Cow_1","Hum_0","Hum_1", "AltBase","AltLength", "BOAEC-D1"), new = c("Fragcount_ref", "Fragcount_alts", "CowMeanscore_ref", "CowMeanscore_alts", "HumMeanscore_ref", "HumMeanscore_alts", "AltBases", "AltLengths", "BAOEC-D1"))
#rm(snps)
#rm(nonsnps)

# Add the multi-allelic variants
multi <- fread(multi_file)
multi_dropcols <- colnames(multi)[str_detect(colnames(multi), "all_")]
multi <- multi[, !..multi_dropcols]
setnames(multi, old=c("count_0", "Cow_0", "Hum_0", "BOAEC-D1"), new=c("Fragcount_ref","CowMeanscore_ref","HumMeanscore_ref", "BAOEC-D1"))

# Correct counts as above, then combine the non-ref fragment counts and mean scores (rounded to 3 dp) into a single column
alt_counts <- c("count_1","count_2","count_3","count_4")
multi[, (c("Fragcount_ref", alt_counts)) := lapply(.SD, `/`, 4), .SDcols = c("Fragcount_ref", alt_counts)]
multi[, Fragcount_alts := do.call(paste, c(.SD, sep = ",")), .SDcols = alt_counts]

cow_scores <- c("Cow_1","Cow_2","Cow_3","Cow_4")
hum_scores <- c("Hum_1","Hum_2","Hum_3","Hum_4")
multi[, (c(cow_scores, hum_scores)) := lapply(.SD, round, 3), .SDcols = c(cow_scores, hum_scores)]
multi[, CowMeanscore_alts := do.call(paste, c(.SD, sep=",")), .SDcols = cow_scores]
multi[, HumMeanscore_alts := do.call(paste, c(.SD, sep=",")), .SDcols = hum_scores]
multi <- multi[,-c(..cow_scores, ..hum_scores, ..alt_counts)]

# Truncate long alleles
multi[, AltBases := vapply(
  strsplit(AltBase, ",", fixed = TRUE),
  \(x) paste(ifelse(nchar(x) >= 50, ">50", x), collapse = ","),
  character(1)
)]
multi[, AltBase := NULL]
multi[,RefBase := if_else(nchar(RefBase) <=50, RefBase, ">50")]
multi <- multi[,VarType := 'multi']

dat_all <- rbindlist(list(dat_all, multi), use.names = T, fill = T)

### Fix so missing genotypes (because filtered out) are consistent NA not mixture of NA and blank
dat_all <- dat_all[,`Bi-83-A` := fifelse(`Bi-83-A` == "", NA, `Bi-83-A`)][,`BAOEC-D1` := fifelse(`BAOEC-D1` == "", NA, `BAOEC-D1`)]

dat_all <- as_tibble(dat_all)

fdrThresh <- 0.05
skew_thresh <- 20

# Remove variants that do not appear in the analysis at all, i.e. those with no result for cow transfections or with high cow skew
dat_all <- dat_all %>% filter(Cow_skew<=skew_thresh, !is.na(Cow_p.adjust))

# Add flag for variants that were not included in all analyses
# IncMain means included in the main analysis -- excluding multiallelic variants that only had two alleles in the SuRE samples and were recoded
# IncHum means included in the human-cow comparison, which excluded variants with failed models in either cell type and with high skew in human cells
dat_all <- dat_all %>% mutate(
  IncMain = fifelse(is.na(Pos), "N", "Y"),
  IncHum = fifelse((is.na(Hum_p.value) | Hum_skew > skew_thresh | is.na(Pos)), "N", "Y")
)

# Add significance flags for each test
dat_all <- dat_all %>%
  mutate(
    Cow_signif = if_else(Cow_p.adjust <= fdrThresh, "Y", "N"),
    Hum_signif = if_else(Hum_p.adjust <= fdrThresh, "Y", "N"),
    Interact_signif = if_else(interact_p.adjust <= fdrThresh, "Y", "N")
)

dat <- dat_all %>% filter(IncMain == "Y") %>% 
  dplyr::select(VarID, Chr, Pos, RefLength, RefBase, AltBases) %>%
  rename_with(tolower)


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
  mutate(Nearest_tss_pos = tss_pos,
         Dist_to_tss     = signed_d,
         Nearest_gene_id = gene_id,
         Nearest_gene_biotype = gene_biotype,
         Nearest_gene_name = gene_name)

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


bw_file <- paste0(data_dir, "merged_RPI2_BOAEC-D1_RPI5.deduplicated_without_bl.bw")

bw <- BigWigFile(bw_file)
## get the maximum score in every window
message("starting PROcap lookup")
start_time <- Sys.time()
summ <- summary(bw, win_gr, type = "max", bpparam = bpparam)   # returns DataFrame with $score
message(sprintf("Run time: %0.2f s",
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
scores <- unlist(summ)$score

dat$PROCap_max <- unlist(summ)$score        # NA means no data → treat as 0 if you like


## Add pangenie frequencies
message("Adding frequency information")
pg_freqsfiles <- list.files(path=data_dir, pattern = "AN_AF.tsv", full.names = T)

names(pg_freqsfiles) <- c("Taurine","Indicine")
col_nms_pg <- c("chr","pos","varid","AF","AN")
pg_freqs <- data.table::rbindlist(lapply(pg_freqsfiles, data.table::fread, col.names = col_nms_pg, na.strings = "."), idcol = "breed")
pg_noid <- pg_freqs[is.na(varid)]
pg_freqs <- pg_freqs[!is.na(varid)]
pg_freqs <- dcast(pg_freqs, formula = varid + chr + pos ~ breed, value.var = c("AF","AN"))
pg_freqs[, chr := as.numeric(chr)]
pg_freqs <- pg_freqs[!is.na(chr)]
pg_freqs <- pg_freqs[, -c("AN_Indicine", "AN_Taurine")]

dat <- left_join(dat, pg_freqs)

# Add the ancestral allele to the SNPs
message("Adding ancetral allele information (SNPs only)")
# read in ancestral alleles
anc_file <- paste0(data_dir, "/snps_ancestral.tsv")
anc <- fread(anc_file, header = F, col.names = c("chr","pos","ancestralAllele"))
# remove any not set to a base (I think usually a gap character)
anc <- anc[str_detect(ancestralAllele, "[ACGT]")]
dat <- data.table(dat)
dat <- merge(dat, anc, all.x=T, by=c("chr","pos"))

# Find which allele matches anc
dat[,whichAnc := fifelse(refbase == ancestralAllele, 0, fifelse(altbases == ancestralAllele, 1, NA))]

# Add column for derived allele frequency
# The AF columns give the frequency of the alt allele.
# So if the ancestral allele is the ref then DAF=AF; if the alt allele is the ref then DAF=1-AF (since this is only biallelics)
dat[,tmpAF_T := as.numeric(AF_Taurine)][,tmpAF_I := as.numeric(AF_Indicine)]
dat[,DAF_Indicine := fifelse(is.na(whichAnc), NA, fifelse(whichAnc == 0, tmpAF_I, (1-tmpAF_I)) )]
dat[,DAF_Taurine := fifelse(is.na(whichAnc), NA, fifelse(whichAnc == 0, tmpAF_T, (1-tmpAF_T)) )]
dat[,c("tmpAF_T", "tmpAF_I", "ancestralAllele","whichAnc"):= NULL]

dat <- dat[,-c("reflength","refbase","altbases")]
setnames(dat, old=c("chr","pos","varid"), new=c("Chr","Pos","VarID"))
setkey(dat, Chr,Pos,VarID)

dat_all <- data.table(dat_all)
setkey(dat_all, Chr,Pos,VarID)
dat_all <- merge(dat_all, dat, all=T)

final_order <- c("Chr","Pos","VarID","IncMain","IncHum", "BAOEC-D1","Bi-83-A", "VarType",
                 "RefBase","AltBases", "RefLength","AltLengths",
                 "Fragcount_ref", "Fragcount_alts","CowMeanscore_ref","CowMeanscore_alts", "HumMeanscore_ref", "HumMeanscore_alts","Cow_skew","Hum_skew",
                 "Cow_estimate","Cow_p.value", "Cow_p.adjust",
                 "Hum_estimate","Hum_p.value", "Hum_p.adjust",
                 "interact_p.value", "interact_p.adjust",
                 "Cow_signif","Hum_signif","Interact_signif"
                 )
setcolorder(dat_all, final_order)

message("Writing main outputs")

# Write output
outfile <- paste0(outdir, "/all_vars.annotated.tsv.gz")
fwrite(dat_all, outfile, quote = F, sep="\t", na="NA")

