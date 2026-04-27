library(tidyverse)
library(data.table)
library(GenomicRanges)
library(VariantAnnotation)

library(ggpubr)
library(patchwork)

# In the Nextflow pipeline the script var_pairs.R identified all pairs of regulatory variants (FDR<=0.05) that appeared together on the same fragment. For pairs with >2 possible haplotypes appearing on at least 10 fragments each, it tested whether a model with both variants was a better fit than one with only the most significant variant.
# However we don't want to test all of these pairs. And also want an idea of what happens to each variant pair that wasn't tested.
# Limit to within 200 bp of each other
# This script takes the paired variants and for each cluster of multiple sites linked together by successive pairs, limits this just to pairs involving a focal variant (rank all paired vars by p-value, work down the list and take pairs with that variant from the wider set of pairs, then take all the variant 2s and remove any additional pairs they appear in from the wider set of pairs. Continue until none remain)
# Then do p-value correction and other processing on the smaller set of pairs

source("analysis_nextflow/data/config.cattle.eddie.R")
winsize <- 200 # find pairs of variants within this bp of each other

# Load all the variants for lookup
vcf_dir <- "~/Documents/SuRE/pangenomes_final/cattle-Tann-pg/bychr/"
var_list <- list()
ranges_list <- list()
for(i in seq_along(chr)) {
  chrid <- chr[i]
  vcffile <- paste0(vcf_dir, chrid, ".cattle-Tann-pg.vcf.gz")
  one_vcf <- readVcf(vcffile)
  one_gr <- rowRanges(one_vcf)
  ranges_list[[i]] <- one_gr
  var_list[[i]] <- data.table(VarID = names(one_gr))
}
var_proc <- rbindlist(var_list)

# Load ASE tested biallelic variants
snpsfile <- "~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/final_stats/all_cattle.snps.tsv.gz"
otherfile <- "~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/final_stats/all_cattle.non_snps.tsv.gz"
fdr_cutoff <- 0.05
snps <- fread(snpsfile, na.strings = "NA")
other <- fread(otherfile, na.strings = "NA")
keepcols <- c("VarID", "Chr", "Pos", "Cow_p.adjust", "Cow_p.value", "Cow_skew", "BOAEC-D1","Bi-83-A")
snps <- snps[,..keepcols][,VarType :="SNP"]
keepcols <- c(keepcols, "VarType")
other <- other[,..keepcols]
tested <- rbindlist(list(snps, other))
tested[ ,`Bi-83-A` := na_if(`Bi-83-A`,"")]
tested[ ,`BOAEC-D1` := na_if(`BOAEC-D1`,"")]
tested <- tested[!is.na(Cow_p.adjust)] # since could be tested in hum not cow

# Get multi for checking
multifile <- "~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/final_stats/all_cattle.multiallelic.tsv.gz"
multi <- fread(multifile, na.strings = "NA")
keepcols <- c("VarID", "Chr", "Pos", "Cow_p.adjust", "Cow_p.value", "Cow_skew", "BOAEC-D1","Bi-83-A")
multi <- multi[,..keepcols][,VarType :="multi"]

# Exclude any which are no-pos or high skew
tested <- tested[,signif := fifelse((Cow_p.adjust <= fdr_cutoff & Cow_skew <= 20 & !is.na(Pos)), T, F)]

sigvars <- tested[signif == T]


# Load in the data from var_pairs.R and find which are within the range
pairdir <- ("~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/variant_pairs/")
find_pairs <- function(chrno) {
  # Find all significant var pairs within 200 bp
  chrid <- chr[chrno]
  chr_gr <- ranges_list[[chrno]]
  chrsig <- sigvars[Chr == chrno]
  
  sig_gr <- chr_gr[names(chr_gr) %in% chrsig$VarID]
  sigranges <- copy(sig_gr)
  # expand ranges of significant variants to winsize either side
  start(sigranges) <- start(sigranges) - winsize
  end(sigranges) <- end(sigranges) + winsize
  # Find overlaps (in both directions along the chr)
  hits <- findOverlaps(sigranges, sig_gr)
  varpairs <- data.table(
    VarID_1 = names(sigranges)[queryHits(hits)],
    VarID_2 = names(sig_gr)[subjectHits(hits)]
  )
  varpairs <- varpairs[VarID_1 != VarID_2]
  
  # Load fragment-based pairs
  pair_analysis <- fread(paste0(pairdir, "all_cattle.", chrno, ".variant_pairs.tsv.gz"), na.strings = "NA")
  # some duplicate rows were accidentally left in the output file, remove these now
  pair_analysis <- pair_analysis[keep != FALSE]
  pair_analysis <- pair_analysis[,-c("keep")]
  
  # Reduce the fragment-based pairs to those within 200 bp and in the filtered signif list
  close_pairs <- merge(varpairs, pair_analysis)
  
  # remove vars with <10 fragments for any haplotype (but not NA)
  few_frags <- close_pairs[(`00` < 10 | `01` < 10 | `10` < 10 | `11` < 10)]
  message(paste0("chr ", chrno, ": removing ", nrow(few_frags), " with <10 fragments for any haplotype"))
  close_pairs <- close_pairs[!(pairID %in% few_frags$pairID)]
  
  # Find the distances
  var1_idx <- match(close_pairs$VarID_1, names(chr_gr))
  var1_gr <- chr_gr[var1_idx]
  var2_idx <- match(close_pairs$VarID_2, names(chr_gr))
  var2_gr <- chr_gr[var2_idx]
  close_pairs$distance <- distance(var1_gr, var2_gr)
  
  return(close_pairs)
}

checked_list <- list()
for(i in seq_along(chr)) {
  chrid <- chr[i]
  pair_results <- find_pairs(i)
  checked_list[[i]] <- pair_results
}

allnear <- rbindlist(checked_list)

# After finding all possible pairs, get the focal variant
paired_vars <- data.table(VarID = unique(c(allnear$VarID_1, allnear$VarID_2)))
paired_vars <- merge(paired_vars, tested[,c("VarID", "Cow_p.value")], by="VarID", all.x=T)
paired_vars <- paired_vars[order(Cow_p.value)]

# Produce new data table for only pairs we will keep
pairs_proc <- copy(allnear)
pairs_proc[, active := TRUE]

# Long form: one row per (VarID, row index)
idx <- rbind(
  pairs_proc[, .(Row = .I, VarID = VarID_1)],
  pairs_proc[, .(Row = .I, VarID = VarID_2)]
)
setkey(idx, VarID)

new_pair_list <- vector("list", length(paired_vars$VarID))
i <- 1L

for (focalvar in paired_vars$VarID) {
  # Find rows involving focalvar via the index
  rows <- idx[VarID == focalvar, Row]
  if (!length(rows)) next
  rows <- rows[pairs_proc[rows, active]]
  if (!length(rows)) next
  focPairs <- pairs_proc[rows]
  focPairs[, tagVar := focalvar]
  new_pair_list[[i]] <- focPairs
  
  donevars <- unique(c(focPairs$VarID_1, focPairs$VarID_2))
  
  # Mark corresponding rows inactive
  pairs_proc[
    VarID_1 %chin% donevars | VarID_2 %chin% donevars,
    active := FALSE
  ]
  
  # Drop those VarIDs from the index
  idx <- idx[!VarID %chin% donevars]
  
  if (!any(pairs_proc$active)) {
    message(paste0("No vars left after ", focalvar, ", done ", i))
    break
  }
  i <- i + 1L
}
new_pair_list <- new_pair_list[seq_len(i - 1L)]
allnear_subset <- rbindlist(new_pair_list)

# Save/load the processed pairs

allnear_subset # reduced set of pairs based on focal var, no filter for haps yet

## Correct p-values and find result
tested_pairs <- as_tibble(allnear_subset) %>% mutate(
  pval_adj = p.adjust(twoall_better, method = 'BH'),
  result = case_when(
    haps < 2 ~ 'fewFrags',
    haps == 2 ~ 'cLD',
    haps > 2 & (failedv | failedb) ~ 'modelFail',
    pval_adj > fdr_cutoff ~ 'ns',
    pval_adj <= fdr_cutoff ~ 'signif'
  )
)
tested_pairs <- tested_pairs %>% dplyr::select(-active)

# Check whether focal var is always var1, if not might need to redo model (original script used corrected p-value and occasionally there were ties; here using raw p-value)
tmp <- tested_pairs %>% filter(VarID_1 != tagVar)
tmp %>% dplyr::count(result)
# in this case all are in cLD so model not fitted anyway and irrelevant which was Var1

# Find outcomes before checking for confounding vars
tested_pairs %>% dplyr::count(result)

focal_vars <- tested_pairs %>% dplyr::count(tagVar, result) %>% 
  pivot_wider(id_cols = tagVar, values_from = n, names_from = result) %>%
  replace_na(list(cLD=0, signif=0, ns=0, fewFrags=0, modelFail=0))
focal_vars_anytest <- focal_vars %>% filter(signif >0 | ns > 0)
focal_vars_anytest %>% filter(signif > 0) %>% nrow()

# Find where any tested pairs have untested vars on the same fragments (>50%) that could be causing the apparent effect (other sig vars on the same fragment should have been accounted for)
# This checks for all tested pairs not just significant
model_pairs <- tested_pairs %>% filter(result == 'signif' | result == 'ns')
mp_details <- merge(data.table(VarID=unique(c(model_pairs$VarID_1, model_pairs$VarID_2))), tested, by='VarID')
model_pairs <- merge(data.table(model_pairs), mp_details[,c("VarID","Chr")], by.x="VarID_1", by.y="VarID")
model_pairs <- model_pairs[,-c( "failedv", "failedb")]

# Go through again to get fragments and find other vars on the fragments
barcode_dir <- "~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/exploring/var_pairs/"
project_script <- "~/Documents/SuRE/scripts/SuRE_SV/analysis/other/plots/project_to_ref.sh"
sigfrag_dir <- "~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/signif_frags/"
subgraph_dir <- "~/Documents/SuRE/pangenomes_final/cattle-Tann-pg/chr_subgraphs/"
setwd(barcode_dir)
pair_regions <- NULL
vars_with_frags <- NULL
species <- 'Cow'
for (chrid in unique(model_pairs$Chr)) {
  message(paste0("doing chr ", chrid))
  fragments <- fread(paste0(sigfrag_dir, chrid, ".all_cattle.signif_frag_chunks.txt.gz"))
  fragments <- fragments[cellSpecies == species]
  chr_pairs <- model_pairs %>% filter(Chr==chrid)
  # For each pair get fragments with both
  chr_frags <- NULL
  for (i in seq_along(chr_pairs$VarID_1)) {
    var1_frags <- fragments[VarID == chr_pairs$VarID_1[i]]
    var2_frags <- fragments[VarID == chr_pairs$VarID_2[i]]
    both_frags <- merge(var1_frags, var2_frags, by=c("input","library","Barcode","cDNA","cellSpecies","Rep"), suffixes = c("_1","_2"))
    both_frags <- dcast(both_frags, formula = input + library + Barcode + VarID_1 + VarAllele_1 + VarID_2 + VarAllele_2 + cellSpecies ~ Rep, value.var="cDNA")
    chr_frags <- bind_rows(chr_frags, both_frags)
  }
  uniqbc <- unique(chr_frags$Barcode)
  bc_list <- paste0(barcode_dir, chrid, ".barcodes.txt")
#  write(uniqbc, file = bc_list)
  
  # Call the script to get the fragment details
  sg <- paste0(subgraph_dir, chrid, ".vg")
  bc_posfile <- paste0(sigfrag_dir, chrid, ".all_cattle.signif_frag_positions.txt.gz")
  bc_cmd <- paste(project_script, "-l", bc_list, "-b", bc_posfile, "-g", sg, "-v vg -p", chrid)
#  system(bc_cmd)
  
  # Load in fragment ref positions and join
  frag_positions <- read_tsv(paste0(barcode_dir, chrid, ".fragment_positions.txt")) %>% 
    dplyr::select(Barcode, RefStart, RefEnd)
  chr_frags <- left_join(chr_frags, frag_positions) %>% mutate(Chr=chrid)
  vars_with_frags <- bind_rows(vars_with_frags, chr_frags)
  
  # For each pair, find max and min
  pair_limits <- chr_frags %>% group_by(VarID_1, VarID_2) %>%
    summarise(minPos = min(RefStart, na.rm = T), maxPos = max(RefEnd, na.rm = T), missing = sum(is.na(RefStart)), allfrags = n()) %>% 
    mutate(Chr = chrid)
  pair_regions <- bind_rows(pair_regions, pair_limits)
}

# Write a regions file to get variants in these regions with bcftools
outdata <- pair_regions %>% ungroup %>% dplyr::select(Chr, minPos, maxPos) %>% 
  arrange(Chr, minPos) %>%
  distinct()
#write_tsv(outdata, file=paste0(barcode_dir, "frag_regions.txt"), quote='none', col_names = FALSE)

# (External) command
#bcftools view ~/Documents/SuRE/pangenomes_final/cattle-Tann-pg/cattle-Tann-pg.vcf.gz -R frag_regions.txt | bcftools query -HH -f '%ID\t%CHROM\t%POS\t%END[\t%GT]\n' > overlap_vars.txt

overlap_vars <- fread(file=paste0(barcode_dir, "overlap_vars.txt"))
bc_with_vars <- unique(vars_with_frags[,c("Barcode", "RefStart", "RefEnd", "Chr")])

# Find which are on each fragment
## Overlaps where POS is inside [RefStart, RefEnd]
pos_hits <- bc_with_vars[
  overlap_vars,
  .(Barcode, RefStart, RefEnd, Chr,
    `#ID`,
    POS = i.POS,         # explicitly keep POS from Df2
    END = i.END,         # explicitly keep END from Df2librar
    `BOAEC_D1:GT`, `Bi_83_A:GT`),
  on = .(Chr = CHROM, RefStart <= POS, RefEnd >= POS),
  nomatch = 0,
  allow.cartesian = TRUE
][,-c("RefStart","RefEnd")]

## Overlaps where END is inside [RefStart, RefEnd]
end_hits <- bc_with_vars[
  overlap_vars,
  .(Barcode, RefStart, RefEnd, Chr,
    `#ID`,
    POS = i.POS,
    END = i.END,`BOAEC_D1:GT`, `Bi_83_A:GT`),
  on = .(Chr = CHROM, RefStart <= END, RefEnd >= END),
  nomatch = 0,
  allow.cartesian = TRUE
][,-c("RefStart","RefEnd")]

## Combine; keep one row per (Df1 row, Df2 row) pair
overlaps <- unique(rbind(pos_hits, end_hits))

# Add to fragments list
frag_overlaps <- merge(vars_with_frags, overlaps, by=c("Barcode","Chr"), all.x=T, allow.cartesian = T)
extra_vars <- frag_overlaps[
  , .N,
  by = .(VarID_1, VarID_2, `#ID`, `POS`)
]
extra_vars <- extra_vars[!((`#ID` == VarID_1) | (`#ID` == VarID_2))]
extra_vars <- merge(extra_vars, pair_regions[,c("VarID_1","VarID_2","allfrags")], by = c("VarID_1","VarID_2"))
extra_vars <- extra_vars[,prop_frags := N / allfrags]
var_lookup <- rbindlist(list(snps[,c("VarID", "Cow_p.value", "Cow_p.adjust", "VarType")], other[,c("VarID", "Cow_p.value", "Cow_p.adjust", "VarType")], multi[,c("VarID", "Cow_p.value", "Cow_p.adjust", "VarType")]))
extra_vars <- merge(extra_vars, var_lookup, by.x="#ID", by.y="VarID", all.x=T)
extra_vars[,cow_signif := fifelse(Cow_p.adjust <= 0.05, "Y","N")]

# Add genotypes (from vcf)
extra_vars <- merge(extra_vars, overlap_vars[,c("#ID","BOAEC_D1:GT","Bi_83_A:GT")])

# Remove any untested because of non-varying genotypes (i.e. not actually vars in the data)
extra_vars <- extra_vars[!(is.na(cow_signif) & `BOAEC_D1:GT` == "1|1" & `Bi_83_A:GT` == "1|1")]

# Look at any not included in original analysis that could have a hidden effect -- untested variants and tested+signficiant multiallelics and multi-now-biallelics (noPos)
tested_check <- tested[, c("VarID")][, is_tested := "Y"]
poss_confounders <- merge(extra_vars, tested_check, all.x = T, by.x="#ID", by.y = "VarID")
poss_confounders <- poss_confounders[is.na(is_tested)]
poss_confounders <- poss_confounders[is.na(cow_signif) | cow_signif == "Y"]

poss_confounders_common <- poss_confounders[prop_frags > 0.5]
poss_conf_vars <- distinct(poss_confounders_common[,c("VarID_1", "VarID_2")])
poss_conf_vars$has_conf = "Y"

# Now link back to original list to find how many
model_pairs_conf <- merge(model_pairs, poss_conf_vars, all.x = T, by=c("VarID_1","VarID_2"))
model_pairs_conf %>% dplyr::count(has_conf, result)
has_poss_conf <- model_pairs_conf %>% dplyr::filter(!is.na(has_conf))
has_poss_conf %>% nrow()

# Check the total results
res_summ <- tested_pairs %>% dplyr::count(result)
rej <- c("cLD","fewFrags","modelFail")
not_modelled <- res_summ %>% filter(result %in% rej) %>% pull(n) %>% sum()

modelled_ok <- model_pairs_conf %>% filter(is.na(has_conf))
modelled_ok %>% dplyr::count(result)
# unique focal vars
modelled_ok %>% filter(result == 'signif') %>% 
  distinct(VarID_1) %>%
  nrow()

# Test for effect of distance
# Test whether the distances in the significant set are greater than in the ns
dists_sig <- modelled_ok %>% filter(result == 'signif') %>%
  pull(distance)
dists_ns <- modelled_ok %>% filter(result == 'ns') %>% pull(distance)
wilcox.test(dists_sig, dists_ns, alternative = 'greater', paired = FALSE)

# update results for plotting
tested_confchecked <- tested_pairs %>%
  mutate(result2 = ifelse(pairID %in% has_poss_conf$pairID, "poss_conf", result))

toplot <- tested_confchecked %>% 
  filter(result2 == 'cLD' | result2 == "ns" | result2 == "signif")
toplot$result2 = factor(toplot$result2, levels = c("cLD","ns", "signif"))

# Plot for Supplementary Fig. 5
p1 <- ggplot(toplot, aes(x=result2, y= distance)) + 
  geom_violin(fill="grey70") + 
  geom_boxplot(width=0.2, fill='grey90') + 
  theme_pubr(base_size = 10, base_family = 'Arial') +
  scale_x_discrete(label = c("complete linkage", "not significant", "significant")) + 
  xlab("") + 
  ylab("Distance between variants (nt)")
p1

figdir <- "~/Documents/SuRE/paper/figures/"
outfile <- paste0(figdir, "/Supplement/Supp_Fig5.svg")
ggsave(outfile, p1,
       width= 160,
       height= 120,
       units = "mm",
       bg= 'transparent')




