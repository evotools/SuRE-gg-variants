library(tidyverse)
library(data.table)
library(patchwork)
library(ggpp)
library(ggpubr)
library(cowplot)
library(ggh4x)
library(scales)

# Get data for the fragment plots in the paper

# Hardcoded file locations, change if needed
frags_dir <- "~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/signif_frags/"
allvars <- fread("~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/sharing_data/reduced/signif_vars.reduced_annot.cow_skew_filtered.tsv.gz") # all *significant* vars
figdir <- "~/Documents/SuRE/paper/figures/fragment_plots/"

axis_labsize <- 5 # pt
annot_size <- 5

# Publication fragment plot
# frag_details now just for this var
frag_plot <- function(onevar, frag_details, var_details, spec='all') {
  variantPos <- var_details$pos
  repcols <- colnames(frag_details)[str_detect(colnames(frag_details), 'Rep')]
  varfrags <- as_tibble(frag_details) %>% select(VarID, VarAllele, library, RefStart, RefEnd, all_of(repcols)) %>%
    dplyr::rename("start"="RefStart", "end"="RefEnd")
  
  # exclude any with fragment position not defined
  varfrags <- varfrags %>% filter(!is.na(start))
  
  varDat <- varfrags %>% pivot_longer(cols=!VarID:end, names_to = "sample", values_to = "value")
  varDat<-varDat %>% mutate(animal = str_replace_all(library, c('BOAEC-D1' = 'B.Taurus', 'Bi-83-A' = 'B.Indicus')),
                            VarAllele = str_replace_all(VarAllele, c('0' = 'Reference', '1' = "Alternative")),
                            VarAllele = fct_relevel(VarAllele, "Reference", "Alternative"))
  varDat <- varDat %>% separate_wider_delim(sample, delim = "_", names=c("Species","Rep"), cols_remove = F)
  varDat <- varDat %>% mutate(Rep = str_replace_all(Rep, c("Rep1" = "Cow 1", "Rep2" = "Cow 2", "Rep3" = "Hum 1", "Rep4" = "Hum 2")))
  
  # keep only species required
  if (spec != 'all') {
    varDat <- varDat %>% filter(Species == spec)
  }
  
  varDat <- varDat %>% mutate(srp_bin=cut_interval(value, n=srp_bins))
  
  # Set beforehand how many bins you want on the x axis and work out which x bin each base position will be in
  densDat<-varDat %>% mutate(id = row_number()) %>%
    reframe(pos = mapply(seq, start, end), .by = -c(start, end)) %>% mutate(density=1)
  densDat.binned <- densDat %>% mutate(xbin=cut(pos, breaks=xbins, labels = F, right=F))
  
  # Reduce so that instead of each fragment having as many rows in each bin as it has bp, to keep only one row from each fragment in each bin
  # **This is not a perfect representation as will mean that if a bin is say 50 bp, a fragment that overlaps only the first 5 bp will contribute as much to the density as a fragment that covers all 50 bp.
  densDat.binned <- densDat.binned %>% distinct(across(-pos), .keep_all = TRUE)
  
  # Do the plot
  p1 <- ggplot(densDat.binned, aes(pos, value) ) +
    geom_bin2d(bins = c(xbins-1,ybins)) +
    scale_fill_viridis_c(trans="log10", guide = guide_colourbar(
      barheight    = unit(10, "mm"),
      barwidth     = unit(2,  "mm")
    )) + 
    theme_pubr(base_size = axis_labsize, base_family = 'Arial', border = TRUE, legend = 'right') +
    geom_vline(xintercept=variantPos,linetype="dashed", colour="grey") +
    ylab("Normalised expression") + xlab(paste0("Chromosome ", chr_to_check, " Position (Mb)")) +
    scale_x_continuous(breaks= extended_breaks(3) , labels = unit_format(unit = "", scale = 1e-6, accuracy = 0.001)) +
    scale_y_continuous(expand = c(0, 0)) + 
    labs(fill='Number of\nfragments') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=1),
          panel.border = element_rect(linewidth = 0.2))
    
    if (spec == 'all') {
      p1 <- p1 + facet_nested(.~Species+VarAllele+Rep, nest_line=element_line(linewidth = 0.1))
    } else {
      p1 <- p1 + facet_nested(.~VarAllele + Rep, nest_line=element_line(linewidth = 0.1))
    }
  
  p1 <- p1 + theme(strip.text.x = element_text(size=axis_labsize),
                   strip.background = element_rect(colour=NA))
  
  #p1
  
  barDat<-varDat %>% group_by(srp_bin, VarAllele, Species,.drop = FALSE) %>% dplyr::count() %>% drop_na() %>%
    ungroup() %>% complete(srp_bin,VarAllele, Species, fill=list(n=0))
  
  barDat<-barDat %>% group_by(VarAllele, Species) %>%
    mutate(var_n = sum(n)) %>%
    mutate(freq = sqrt(n / var_n)) %>% ungroup() %>%  
    mutate(freq= ifelse(VarAllele =="Alternative",freq,freq*-1))
  
  p2 <- ggplot(barDat) + 
    geom_bar(aes(y=srp_bin, x=freq, fill=as.factor(VarAllele)),stat="identity",position="identity", width=1) +
    scale_x_continuous(          # <- control axis here
      limits = c(-1, 1),
      breaks = c(-1, 0, 1)       # only these labels will be shown
    ) +
    scale_fill_brewer(palette = "Accent", name = "Allele") +
    theme_pubr(base_size = axis_labsize, base_family = 'Arial', border = TRUE) +
    xlab("Square root\nproportion") +
    theme(panel.grid.major.x = element_line(colour = 'darkgrey', linewidth = 0.1),
          panel.grid.minor.x = element_line(colour = 'darkgrey', linewidth = 0.1),
          axis.title.y = element_blank(),
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks=element_blank(),
          legend.text = element_text(size = annot_size),
          legend.title = element_text(size = annot_size),
          legend.key.width = unit(2, "mm"),
          legend.key.height = unit(2, "mm"),
          panel.border = element_rect(linewidth = 0.2)
    )
  if (spec == 'all') {
    p2 <- p2 + facet_wrap(~Species) + theme(strip.text.x = element_text(size=axis_labsize),
                                            strip.background = element_rect(colour=NA))
  }
  return(list(p1=p1, p2=p2))
}

# Figure 2C and 5C
# SNV near CLSTN2, with and without human cell data

clstn2_id <- ">1227050>1227053"
thisvar_details <- allvars[varid==clstn2_id]
chr_to_check <- 1
# Get the fragment details for this variant
chr_bc <- fread(paste0(frags_dir, chr_to_check, ".all_cattle.signif_frag_scores.txt.gz"))
var_bc <- chr_bc[VarID == clstn2_id][,chr := chr_to_check]


# Get the fragment positions
# These have been extracted using project_to_ref.sh
refpos <- fread("~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/signif_vars/topN/toref/1.fragment_positions.txt")
refpos <- refpos[,-c("Start","End","RefStrand")]
frag_details_refpos <- merge(var_bc, refpos, by = "Barcode", all.x = T)
# check for missing -- none
frag_details_refpos[is.na(RefStart)]
# save the details
save(frag_details_refpos, thisvar_details, clstn2_id, chr_to_check, file = paste0(figdir, "CLSTN2/CLSTN2.fragdetails.Rdata"))

# or else just load the pre-prepared info
#load(paste0(figdir, "CLSTN2/CLSTN2.fragdetails.Rdata"))

plot_props <- c(4, 0.9) # for output, relative proportions of the fragment plot and srp plot)
spacer_width <- 0.2 # relative space between plots
#plot_height <- 60 # mm
#plot1_width <- 100 # width of output for single set of plots

# Bin sizes can change between variants
ybins <- 70 # srp have one block per 2 ybins
srp_bins <- ybins/2
xbins <- 70

# Make the plot for cow only
#clstn_cow <- frag_plot(clstn2_id, frag_details_refpos, thisvar_details, 'Cow')
clstncow_plots <- frag_plot(clstn2_id, frag_details_refpos, thisvar_details, 'Cow')
cp1 <- clstncow_plots[['p1']]
cp2 <- clstncow_plots[['p2']]
cow_pw <- cp1 + cp2 + plot_layout(widths=plot_props, guides = "collect") & theme(legend.position = 'right')
print(cow_pw)

# Save image
outfile <- paste0(figdir, "/CLSTN2/CLSTN2_cow.svg")
ggsave(outfile, cow_pw,
       width = 95,
       height = 50,
       units = "mm",
       bg= 'transparent')

# Save for putting in multipanel figure
save(cp1, cp2, cow_pw, file="~/Documents/SuRE/paper/figures/Fig2/fragplot.Rdata")

# Plot for human
clstnhum_plots <- frag_plot(clstn2_id, frag_details_refpos, thisvar_details, 'Hum')
hp1 <- clstnhum_plots[['p1']]
hp2 <- clstnhum_plots[['p2']]

# Work out how much width should increase to add second plots and spacer
# orig_width <- plot_props[1] + plot_props[2]
# double_width <- orig_width*2 + spacer_width
# plot2_width <- double_width/orig_width * plot1_width
# # This isn't working, I think because of the legend
# plot2_width <- 1.95*plot1_width

# Actually deliberately use different sizing -- showing something different, doesn't need to be identical. First plot should be larger as introducing the plot layout

clstn_both <- cp1 + cp2 + plot_spacer() + hp1 + hp2 + plot_layout(widths=c(plot_props, spacer_width, plot_props), guides = "collect", nrow = 1, ncol = 5) & theme(legend.position = 'right')
clstn_both
outfile <- paste0(figdir, "/CLSTN2/CLSTN2_both.svg")
ggsave(outfile, clstn_both,
       width = 180,
       height = 50,
       units = "mm",
       bg= 'transparent')




# Figure 6D -- xsQTL
# SNV near BMP8B but closer to a gene not annotated but lifted over from human. So just call it xstql_eg
xs_qtl <- ">43578898>43578901"
xsqtl_details <- allvars[varid==xs_qtl]
chr_to_check <- 3
# Get the fragment details for this variant
chr_bc <- fread(paste0(frags_dir, chr_to_check, ".all_cattle.signif_frag_scores.txt.gz"))
var_bc <- chr_bc[VarID == xs_qtl][,chr := chr_to_check]

# Get the fragment positions
# These have been extracted using project_to_ref.sh
xsqtl_refpos <- fread("~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/signif_vars/topN/toref/3.fragment_positions.txt")
xsqtl_refpos <- xsqtl_refpos[,-c("Start","End","RefStrand")]
frag_details_xsqtl_refpos <- merge(var_bc, xsqtl_refpos, by = "Barcode", all.x = T)
# check for missing -- none
frag_details_xsqtl_refpos[is.na(RefStart)]
# save the details
save(frag_details_xsqtl_refpos, xsqtl_details, xs_qtl, chr_to_check, file = paste0(figdir, "xsQTL/xsQTL.43578898_43578901.fragdetails.Rdata"))
load(paste0(figdir, "xsQTL/xsQTL.43578898_43578901.fragdetails.Rdata"))


# Make the plots
ybins <- 70 # srp have one block per 2 ybins
srp_bins <- ybins/2
xbins <- 70
plot_props <- c(4, 0.8)

# Cow
xsqtl_cowplots <- frag_plot(xs_qtl, frag_details_xsqtl_refpos, xsqtl_details, 'Cow')
xs_cp1 <- xsqtl_cowplots[['p1']]
xs_cp2 <- xsqtl_cowplots[['p2']]
# Hum
xsqtl_humplots <- frag_plot(xs_qtl, frag_details_xsqtl_refpos, xsqtl_details, 'Hum')
xs_hp1 <- xsqtl_humplots[['p1']]
xs_hp2 <- xsqtl_humplots[['p2']]

xsqtl_both <- xs_cp1 + xs_cp2 + plot_spacer() + xs_hp1 + xs_hp2 + 
  plot_layout(widths=c(plot_props, spacer_width, plot_props), guides = "collect", nrow = 1, ncol = 5) &
  theme(legend.position = 'right')

# Both together
xsqtl_bothplot <- frag_plot(xs_qtl, frag_details_xsqtl_refpos, xsqtl_details, 'all')
xsqtl_both2 <- xsqtl_bothplot[['p1']] + plot_spacer() + xsqtl_bothplot[['p2']] +
  plot_layout(widths=c(plot_props[1], spacer_width, plot_props[2]), guides = "collect", nrow = 1, ncol = 3) &
  theme(legend.position = 'right')
xsqtl_both2
xsqtl_both <- xsqtl_both2


# Or better to show both at the same time?
# xsqtl_bothplots <- frag_plot(xs_qtl, frag_details_xsqtl_refpos, xsqtl_details, 'all')
# xs_p1 <- xsqtl_bothplots[['p1']]
# xs_p2 <- xsqtl_bothplots[['p2']]
# xsqtl_both <- xs_p1 + xs_p2 + plot_layout(widths=plot_props, guides = "collect") & theme(legend.position = 'right')
xsqtl_both


outfile <- paste0(figdir, "/xsQTL/xsQTL_both.svg")
ggsave(outfile, xsqtl_both,
       width = 180,
       height = 55,
   #    height = 50,
       units = "mm",
       bg= 'transparent')




### Fig 2e
# Multi varplot
figdir <- "~/Documents/SuRE/paper/figures/fragment_plots/multi/"
frags_dir_multi <- "~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/signif_frags_multi/"
multi_ex <- ">48800609>48800613"
multi_details <- allvars[varid==multi_ex]
chr_to_check <- 5
# Get the fragment details for this variant
chr_bc <- fread(paste0(frags_dir_multi, chr_to_check, ".all_cattle.multi.signif_frag_scores.txt.gz"))
var_bc <- chr_bc[VarID == multi_ex][,chr := chr_to_check]

multi_refpos <- fread("~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/signif_vars/topN_multi/byref/5.fragment_positions.txt")
multi_refpos <- multi_refpos[,-c("Start","End","RefStrand")]
frag_details_multi_refpos <- merge(var_bc, multi_refpos, by = "Barcode", all.x = T)
# check for missing
frag_details_multi_refpos[is.na(RefStart)]
# save the details
save(frag_details_multi_refpos, multi_details, multi_ex, chr_to_check, file = paste0(figdir, "multi_data.Rdata"))
load(paste0(figdir, "multi_data.Rdata"))

frag_plot_multi <- function(onevar, frag_details, var_details, spec='all') {
  variantPos <- var_details$pos
  repcols <- colnames(frag_details)[str_detect(colnames(frag_details), 'Rep')]
  varfrags <- as_tibble(frag_details) %>% select(VarID, VarAllele, library, RefStart, RefEnd, all_of(repcols)) %>%
    dplyr::rename("start"="RefStart", "end"="RefEnd")
  
  # exclude any with fragment position not defined
  varfrags <- varfrags %>% filter(!is.na(start))
  
  varDat <- varfrags %>% pivot_longer(cols=!VarID:end, names_to = "sample", values_to = "value")
  
  alleles <- unique(varDat$VarAllele)
  allele_names <- c('0' = 'Reference', '1' = "Alternative 1", "2" = "Alternative 2", "3" = "Alternative 3", "4" = "Alternative 4")
  allele_names <- allele_names[names(allele_names) %in% alleles]
  
  varDat<-varDat %>% mutate(animal = str_replace_all(library, c('BOAEC-D1' = 'B.Taurus', 'Bi-83-A' = 'B.Indicus')),
                            VarAllele = str_replace_all(VarAllele, allele_names),
                            VarAllele = fct_relevel(VarAllele, allele_names))
  varDat <- varDat %>% separate_wider_delim(sample, delim = "_", names=c("Species","Rep"), cols_remove = F)
  varDat <- varDat %>% mutate(Rep = str_replace_all(Rep, c("Rep1" = "Cow 1", "Rep2" = "Cow 2", "Rep3" = "Hum 1", "Rep4" = "Hum 2")))
  
  # keep only species required
  if (spec != 'all') {
    varDat <- varDat %>% filter(Species == spec)
  }
  
  
  # Set beforehand how many bins you want on the x axis and work out which x bin each base position will be in
  densDat<-varDat %>% mutate(id = row_number()) %>%
    reframe(pos = mapply(seq, start, end), .by = -c(start, end)) %>% mutate(density=1)
  densDat.binned <- densDat %>% mutate(xbin=cut(pos, breaks=xbins, labels = F, right=F))
  
  # Reduce so that instead of each fragment having as many rows in each bin as it has bp, to keep only one row from each fragment in each bin
  # **This is not a perfect representation as will mean that if a bin is say 50 bp, a fragment that overlaps only the first 5 bp will contribute as much to the density as a fragment that covers all 50 bp.
  densDat.binned <- densDat.binned %>% distinct(across(-pos), .keep_all = TRUE)
  
  # Do the plot
  p1 <- ggplot(densDat.binned, aes(pos, value) ) +
    geom_bin2d(bins = c(xbins-1,ybins)) +
    scale_fill_viridis_c(trans="log10", guide = guide_colourbar(
      barwidth    = unit(10, "mm"),
      barheight     = unit(2,  "mm")
    )) + 
    theme_pubr(base_size = axis_labsize, base_family = 'Arial', border = TRUE, legend = 'top') +
    geom_vline(xintercept=variantPos,linetype="dashed", colour="grey") +
    ylab("Normalised expression") + 
    xlab(paste0("Chromosome ", chr_to_check, " Position (Mb)")) +
    scale_x_continuous(breaks= extended_breaks(3) , labels = unit_format(unit = "", scale = 1e-6, accuracy = 0.001)) +
    scale_y_continuous(expand = c(0, 0)) + 
    labs(fill='Number of\nfragments') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=1),
          panel.border = element_rect(linewidth = 0.2))
  
  if (spec == 'all') {
    p1 <- p1 + facet_nested(.~Species+VarAllele+Rep, nest_line=element_line(linewidth = 0.1))
  } else {
    p1 <- p1 + facet_nested(.~VarAllele + Rep, nest_line=element_line(linewidth = 0.1))
  }
  
  p1 <- p1 + theme(strip.text.x = element_text(size=axis_labsize),
                   strip.background = element_rect(colour=NA))
  
  p1
  return(p1)
}

ybins <- 70 
xbins <- 70

multi_fragplot <- frag_plot_multi(multi_ex, frag_details_multi_refpos, multi_details, spec='Cow')
multi_fragplot

save(multi_fragplot, file="~/Documents/SuRE/paper/figures/fragment_plots/multi/multi_fragplot.Rdata")
