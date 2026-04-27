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
figdir <- "~/Documents/SuRE/paper/figures/fragment_plots/indel"

axis_labsize <- 5 # pt
annot_size <- 5


##### Node plot #####
indel_id <- ">7092778>7092780"
indel_details <- allvars[varid==indel_id]
chr_to_check <- indel_details$chr

# Get the fragment details for this variant
chr_bc <- fread(paste0(frags_dir, chr_to_check, ".all_cattle.signif_frag_scores.txt.gz"))
var_bc <- chr_bc[VarID == indel_id][,chr := chr_to_check]
nodelen_dir <- '~/Documents/SuRE/pangenomes_final/cattle-Tann-pg/chr_subgraphs/node_lengths/'
chr_nodelens <- read_tsv(file=paste0(nodelen_dir, chr_to_check, ".node_lengths.txt"), col_names = c("nodeID",'nodelength'))

path_colnames <- c("Barcode", "Start","End","Haplotype","library","fullPath")
pathfiles <- list.files("~/Documents/SuRE/analysis/cattle/cattle-Tann-pg/signif_vars/topN/paths/", pattern = paste0("\\.",chr_to_check, "\\.frag_paths.txt"), full.names = T)
path_list <- lapply(pathfiles, function(file) {
  fread(file, col.names = path_colnames, header = FALSE)
})
fragpaths <- rbindlist(path_list)
fragpaths <- fragpaths[,library := NULL]
rm(path_list)

frag_details_path <- merge(var_bc, fragpaths, by = "Barcode", all.x = T)

# Read paths could be missing because the read paths had multiple matches in the haplotype (rare) or because the read path was not exactly matched (still rare but more common)
# Paths (or rather node lists) may also contain a ? -- this is where the read paths matched both haplotypes, and the path between them differed between haplotypes. Nodes found in both haplotypes are given, with a ? in the path where it could be one of several possibilities. More common than above. But by definition there must be some non-? nodes (otherwise it would have been possible to distinguish them)
pathmiss <- frag_details_path %>% select(VarID, VarAllele, fullPath) %>%
  mutate(match = case_when(
    is.na(fullPath) ~ 'multi',
    fullPath == "NotMatched" ~ 'nomatch',
    str_detect(fullPath, "\\?") ~ 'ambig',
    .default = "ok"
  ))


# save
save(indel_id, indel_details, frag_details_path, pathmiss, file = paste0(figdir, "/indel_fragdetails.Rdata"))

# reload
load(paste0(figdir, "/indel_fragdetails.Rdata"))

# set parameters
nodeplot_params <- list()
nodeplot_params$flanksize <- 6  # How many nodes either side of the variant to plot
nodeplot_params$padding <- 5  # small gap between the max yvalue and edge of the plot to make sure it's drawn properly
nodeplot_params$ybinsize <- 4
nodeplot_params$gap <- 10 # x2 to get the height of the gap in the middle between the upper and lower plot where the nodes are drawn
nodeplot_params$fig.height <- nodeplot_params$gap - 5 # size of the node rectangles within the gap
nodeplot_params$twocols <- c('grey30', 'grey80') # alternating colours for the alternating nodes
nodeplot_params$ylabel.xspace <- 0.2  # distance of y label from x axis as a fraction of the total xwidth

flanksizes <- c(10,3)

list2env(nodeplot_params, envir = environment())

# Function to extract the y start and end from the automatically named bins
bin_limits <- function(binlabel) {
  binlabel <- str_replace_all(binlabel, "[\\[\\]\\(\\)]", "")
  limits <- as.numeric(str_split(binlabel,",")[[1]])
  
  return(list(ystart=limits[1], yend=limits[2]))
}

# Process data for the plot
spec <- 'Cow'
repcols <- colnames(frag_details_path)[str_detect(colnames(frag_details_path), 'Rep')]
fragments <- as_tibble(frag_details_path) %>% select(chr, library, VarAllele, fullPath, all_of(repcols))
fragments <- fragments %>% filter(!is.na(fullPath)) %>%
  filter(fullPath != "NotMatched")
fragments <- fragments %>% pivot_longer(!chr:fullPath, names_to='sample', values_to = 'value')
fragments <- fragments %>% separate_wider_delim(sample, delim = "_", names=c("Species","Rep"), cols_remove = F)
fragments <- fragments %>% mutate(Rep = str_replace_all(Rep, c("Rep1" = "Cow 1", "Rep2" = "Cow 2", "Rep3" = "Hum 1", "Rep4" = "Hum 2")))

# keep only species required
if (spec != 'all') {
  fragments <- fragments %>% filter(Species == spec)
}
fragments <- fragments %>% mutate(VarAllele = str_replace_all(VarAllele, c('0'='Reference', '1'= 'Alternative')), VarAllele = fct_relevel(VarAllele, "Reference", "Alternative"))

varnodes <- as.numeric(strsplit(indel_id, "<|>")[[1]][2:3])

# Split into nodes
fragments.nodes <- fragments %>% separate_longer_delim(fullPath, delim=",") %>%
  filter(fullPath != '?') %>%
  mutate(node=as.numeric(fullPath), .keep='unused')

plotnode1 <- varnodes[1] - flanksize
plotnode2 <- varnodes[2] + flanksize
frags.plot <- fragments.nodes %>% filter(node >= plotnode1, node <=plotnode2)

densDat.binned <- frags.plot %>% mutate(ybin= if_else(VarAllele=="Reference", cut_interval(value, length=ybinsize, right=F), cut_interval(-value, length=ybinsize)))
densDat.binned <- densDat.binned %>% mutate (new_cols = lapply(ybin, bin_limits)) %>%
  unnest_wider(new_cols)
densDat.tiled <- densDat.binned %>% group_by(node,ybin,Rep,VarAllele, ystart, yend) %>% summarise(n=n()) %>% ungroup()
densDat.tiled <- densDat.tiled %>% left_join(chr_nodelens, by=c("node"="nodeID"))

# Get node positions and colours
xwidths <- densDat.tiled %>% select(node, nodelength) %>% distinct() %>% arrange(node)
xwidths <- xwidths %>% mutate(xstart = cumsum(lag(nodelength, default=0)) + 1, xend=xstart+nodelength)

xwidths <- xwidths %>% mutate(colour=rep(1:2, length.out=n()))
xwidths1 <- xwidths %>% filter(colour==1)
xwidths2 <- xwidths %>% filter(colour==2)

densDat.tiled <- densDat.tiled %>% left_join(xwidths)
densDat.tiled <- densDat.tiled %>% mutate(plotystart = if_else(VarAllele=='Reference', ystart+gap, ystart - gap), plotyend = if_else(VarAllele=='Reference', yend+gap, yend - gap))

# Get correct y axis tickmarks taking the middle figure into account
# Both upper and lower plot go to the same max, whichever is higher
ymax.val <- max(c(max(densDat.tiled$yend), abs(min(densDat.tiled$ystart))))
ymax <- ymax.val + gap + padding
if (ymax.val < 0.5) {
  ymax <- 0.5 + gap + padding
}
ymin <- -ymax

# Set y-axis labels
labelpos.all <- c(rev(seq(from=-gap, to=-1000-gap, by=-50)), seq(from=gap, to=1000+gap, by=50))
labels.all <- c(seq(from=1000, to=0, by=-50), seq(from=0, to=1000, by=50))
labelpos <- labelpos.all[which(labelpos.all >= ymin & labelpos.all <= ymax)]
labels <- labels.all[which(labelpos.all >= ymin & labelpos.all <= ymax)]

# Add y axis label for each allele
# This will be added to the first replicate of Cow only - edit if using all or human only
datawidth <- max(xwidths$xend)
ylabel.xpos <- 0-(datawidth * ylabel.xspace) - datawidth * 0.05
ylabel <- tibble(ypos = c(ymax/2 + gap, 0-ymax/2 - gap), xpos = ylabel.xpos, Rep="Cow 1", text=c("Normalised expression\n(reference)", "Normalised expression\n(alternative)"))

# Find y bins that have non-0 for some x values but not others and add dashed lines
# if sample has ystart and y end values then get a row for all possible xstart-xend
to.add <- densDat.tiled %>% group_by(plotystart, plotyend, Rep) %>% mutate(minstart = min(xstart), maxend=max(xend)) %>% ungroup %>%
  complete(nesting(plotystart, plotyend, Rep, minstart, maxend), nesting(xstart,xend)) %>% 
  filter(is.na(ybin), xstart > minstart, xend < maxend) %>% select(plotystart:xend) %>%
  mutate(yval=(plotystart+plotyend)/2)
tmp <- densDat.tiled %>% group_by(plotystart, plotyend) %>% mutate(minstart = min(xstart), maxend=max(xend)) %>% ungroup %>% 
  complete(nesting(plotystart, plotyend, Rep, minstart, maxend), nesting(xstart,xend)) %>% 
  mutate(yval=(plotystart+plotyend)/2)


# where to indicate variant (position comes from the ends of the presumed flanking nodes as given in the variant ID)
ylimits <- densDat.tiled %>% group_by(Rep) %>% summarise(ymax=max(plotyend)+3, ymin =min(plotystart)-3)
varlabel <- xwidths %>% filter (node == varnodes[1] | node == varnodes[2]) %>% summarise(xmin = min(xend), xmax=max(xstart))
varlabel.df <- tibble(ylimits, xmin=varlabel$xmin, xmax=varlabel$xmax)
varlabel.df <- varlabel.df %>% pivot_longer(c(xmin, xmax), names_to = "x_type", values_to = "xvalues") %>%
  pivot_longer(c(ymin, ymax), names_to = "y_type", values_to = "yvalues")
tmp1 <- varlabel.df %>% filter(x_type=='xmax') %>% arrange(y_type)
tmp2 <- varlabel.df %>% filter(x_type=='xmin') %>% arrange(desc(y_type))
varlabel.df <- bind_rows(tmp1, tmp2, tmp1[1:2,])


fragplot <- 
  ggplot( ) +
  geom_rect(densDat.tiled , mapping= aes(xmin=xstart, xmax=xend, ymin=plotystart, ymax=plotyend, fill=n)) + 
  scale_fill_viridis_c(trans="log10", guide = guide_colourbar(
    barheight    = unit(2, "mm"),
    barwidth     = unit(10,  "mm")
  ),
  name = "Number of\nfragments",
  na.value = 'transparent') + 
  theme_pubr(base_size = axis_labsize, base_family = 'Arial', border = TRUE, legend = 'bottom') +
  geom_hline(yintercept = gap) + 
  geom_hline(yintercept = -gap) +
  scale_y_continuous(breaks=labelpos, labels=labels, limits=c(ymin,ymax)) +
  xlab('Cumulative node length (nt)') +
  ylab("") +
  facet_wrap(~Rep) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        panel.border = element_rect(linewidth = 0.2),
        strip.text.x = element_text(size=axis_labsize),
        strip.background = element_rect(colour=NA),
   #     plot.margin = margin(2,2,2,8, unit = 'mm'),
        axis
  ) +
  coord_cartesian(clip="off", xlim=c(0,datawidth*1.02))

fragplot_annot <- fragplot + 
  annotate('rect', xmin=xwidths1$xstart, xmax=xwidths1$xend, ymin=-fig.height, ymax=fig.height, fill=twocols[1]) + 
  annotate('rect', xmin=xwidths2$xstart, xmax=xwidths2$xend, ymin=-fig.height, ymax=fig.height, fill=twocols[2]) +
  geom_segment(data=to.add, mapping=aes(x=xstart, y=yval, xend = xend, yend = yval), colour='grey80', linetype='solid', linewidth =0.2) +
  geom_text(data = ylabel, mapping=aes(x=xpos, y=ypos, label=text, angle=90, vjust=0)) +
  geom_path(data=varlabel.df, mapping=aes(x=xvalues, y=yvalues), colour='red', linetype='dashed', linewidth=0.2)

fragplot_annot



# output plots
save(fragplot_annot, file=paste0(figdir, "/indel_nodeplot.Rdata"))

outfile <- paste0(figdir, "/indel_nodeplot.svg")
ggsave(outfile, fragplot.annot,
       width = 90,
       height = 60,
       units = "mm",
       bg= 'transparent')
