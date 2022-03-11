#!/usr/bin/env Rscript

# comp_mut_tree.R

# Description
# Plot a newick tree file with ggtree
# ggtree tutorial - https://garynapier.github.io/ggtree/
# Install ggtree - https://bioconductor.org/packages/release/bioc/html/ggtree.html (code: BiocManager::install("ggtree"))

# Arguments to script:

# Input:

# Main steps:

# Output:

# RUN:
# Rscript r_scripts/<file_name>.R
# Rscript r_scripts/<file_name>.R

# Setup ----

library(optparse)
library(phytools)
library(ggtree)
library(scales)
library(gplots)
library(castor)
library(ggnewscale)

source("https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Functions.R")

# # Arguments ----
# 
# option_list = list(
#   # EXAMPLE:
#   # make_option(c("-t", "--template_file_name"), type="character", default=NULL,
#   #             help="input template xml file", metavar="character"),
#   
#   make_option(c("-", "--"), type="", default=NULL, 
#               help="", metavar="")
# ); 
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# print("ARGUMENTS:")
# print(opt)
# print("---")
# print(str(opt))
# 
# # Files ----
# 
# template_file <- opt$template_file_name
# 
# print("FILES:")
# print(c(": ", ))

setwd("~/Documents/comp_mut/")

# Setup
circular <- T
if(circular){
  layout <- "circular"
}else{
  layout <- "rectangular"
}

project_code <- 'isoniazid'
project_code_col_name <- 'drug'

# Files ----

tree_file <- "newick/isoniazid.filt.val.gt.g.snps.fa.treefile"
# resistance_mutations_file <- "results/potential_res_mut_samps.csv"
metadata_file <- "results/potential_res_mut_samps.csv"

# Read in data ----

tree_all_samps <- read.newick(tree_file)
# res_mut <- read.csv(resistance_mutations_file, header = T)
metadata <- read.csv(metadata_file, header = T)

# Wrangle data ----

# Midpoint root the tree
tree_all_samps <- phytools::midpoint.root(tree_all_samps)
n_samps <- length(tree_all_samps$tip.label)

# Metadata
# Subset metadata to just the project code rows
metadata <- metadata[metadata[, project_code_col_name] == project_code, ]
row.names(metadata) <- metadata$wgs_id
# Lineages - remove 'lineage' and convert to factor
metadata$main_lineage <- factor(gsub('lineage', '', metadata$main_lineage))
metadata$sublin <- factor(gsub('lineage', '', metadata$sublin))


lin_data <- metadata[, "main_lineage", drop = F]
dr_status_data <- metadata[, "drtype", drop = F]
dr_data <- metadata[,'dst', drop = F]
lin_data <- data.frame(apply(lin_data, 2, as.factor))
dr_status_data <- data.frame(apply(dr_status_data, 2, as.factor))
dr_data <- data.frame(apply(dst_data, 2, as.factor))
# row.names(lin_data) <- metadata$wgs_id
# row.names(dr_status_data) <- metadata$wgs_id

# Colours
alpha <- 0.9
lin_colours <- rainbow(length(unique(metadata$main_lineage)), alpha = alpha)
dr_df <- data.frame(drtype = c("Sensitive", "Pre-MDR-TB", "MDR-TB", "Pre-XDR-TB", "XDR-TB", "Other"), 
                    col = c("green1", "yellow2", "orange1", "red1", "black", "grey"))
dr_status_colours <- scales::alpha(gplots::col2hex(dr_df$col), alpha = alpha)
names(lin_colours) <- c(sort(unique(metadata$main_lineage)))
names(dr_status_colours) <- dr_df$drtype

# Set up ggtree parameters ----

width <- 0.05
font_sz <- 3
line_sz <- 0.25
angle <- 30

y_lim <- c(-10, n_samps + (n_samps * 0.1))
legend_spec <- theme(legend.title = element_text(size = 9),
                     legend.text = element_text(size = 7),
                     legend.key.size = unit(0.3, "cm"))

max_dist <- castor::get_tree_span(tree_all_samps, as_edge_count=FALSE)$max_distance
# min_dist <- castor::get_tree_span(tree_all_samps, as_edge_count=FALSE)$min_distance
# offset <- max_dist

# Make basic tree ----

ggtree_all_samps <- ggtree(tree_all_samps, 
                           size = line_sz)+
  coord_cartesian(ylim = y_lim)

# Add lineage data 
lin_hm <- gheatmap(ggtree_all_samps, lin_data, 
                   width = width, 
                   offset = 0, 
                   colnames_position = "top",
                   colnames_angle = angle, 
                   colnames_offset_y = 1,
                   hjust = 0,
                   font.size = font_sz) +
  # geom_vline(aes(xintercept = max_dist), col = "red")+
  # Add the custom colours defined above
  scale_fill_manual(values = lin_colours, breaks = names(lin_colours) ) +
  # Define the legend title
  labs(fill = "Lineage")

ggtree_data <- ggplot_build(lin_hm)
# pos <- unique(x$data[[3]]$x)
wd <- unique(ggtree_data$data[[3]]$xmax)-unique(ggtree_data$data[[3]]$xmin)

lin_hm <- lin_hm + ggnewscale::new_scale_fill() 

# Add DR status
dr_status_hm <- gheatmap(lin_hm, dr_status_data,
                         width = width,
                         offset = wd,
                         color = NULL,
                         colnames_position = "top",
                         colnames_angle = angle, 
                         colnames_offset_y = 1,
                         hjust = 0,
                         font.size = font_sz) +
  # geom_vline(aes(xintercept = pos), col = "red")+
  scale_fill_manual(values = dr_status_colours, breaks = names(dr_status_colours) )+
  labs(fill = "DR\nstatus")

dr_status_hm

# Add DST data


dr_status_hm <- dr_status_hm + ggnewscale::new_scale_fill() 

gheatmap(dr_status_hm, dr_data,
         # Increase offset
         offset = wd*2,
         width = width,
         # Change color to black
         # color = NULL,
         color="black",
         low="white", 
         high="black", 
         colnames_position = "top",
         colnames_angle = angle, 
         colnames_offset_y = 1,
         hjust = 0,
         font.size = 2.5) +
  # Define colours
  scale_fill_manual(values=c("white", "black"), labels = c("Sensitive", "Resistant", "NA"), na.value = "grey")+
  labs(fill = "Drug\nresistance")











