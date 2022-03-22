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
library(ggplot2)
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
outfile <- paste0("results/", project_code, "_tree.png")

# Read in data ----

tree_all_samps <- read.newick(tree_file)
# res_mut <- read.csv(resistance_mutations_file, header = T)
metadata <- read.csv(metadata_file, header = T)

# Wrangle data ----

drug_abv_table <- data.frame(drug_full_nm = c("isoniazid", "rifampicin"),
                             drug_abv = c("INH", "RIF"))

drug_abv <- drug_abv_table[drug_abv_table["drug_full_nm"] == project_code, "drug_abv"]

# Midpoint root the tree
tree_all_samps <- phytools::midpoint.root(tree_all_samps)
n_samps <- length(tree_all_samps$tip.label)

# Metadata
# Subset metadata to just the project code rows
metadata <- metadata[metadata[, project_code_col_name] == project_code, ]
# Make ID rownames of metadata so heatmap strips work in ggtree
row.names(metadata) <- metadata$wgs_id
# Lineages - remove 'lineage' and convert to factor
metadata$main_lineage <- factor(gsub('lineage', '', metadata$main_lineage))
metadata$sublin <- factor(gsub('lineage', '', metadata$sublin))
# Convert DST to factor
metadata$dst <- as.factor(metadata$dst)

# Get data for separate heatmap strips
lin_data <- metadata[, "main_lineage", drop = F]
dr_status_data <- metadata[, "drtype", drop = F]
dr_data <- metadata[,'dst', drop = F]

# Change col headers to match legends 
lin_data_lab <- "Lineage"
dr_status_lab <- "DR status"
dr_data_lab <- paste(drug_abv, "DST")

colnames(lin_data) <- lin_data_lab
colnames(dr_status_data) <- dr_status_lab
colnames(dr_data) <- dr_data_lab

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

# Make tree ----

ggtree_all_samps <- ggtree(tree_all_samps, 
                           size = line_sz)
  # coord_cartesian(ylim = y_lim)

# Add lineage data 
lin_hm <- gheatmap(ggtree_all_samps, lin_data, 
                   width = width, 
                   offset = 0, 
                   colnames_position = "top",
                   colnames_angle = angle, 
                   colnames_offset_y = 1,
                   hjust = 0,
                   font.size = font_sz) +
  # Add the custom colours defined above
  scale_fill_manual(values = lin_colours, breaks = names(lin_colours) ) +
  # Define the legend title
  labs(fill = lin_data_lab)

# Pull the width of the strip from the plot just created in order to set offset for next strips
ggtree_data <- ggplot2::ggplot_build(lin_hm)
os <- unique(ggtree_data$data[[3]]$xmax)-unique(ggtree_data$data[[3]]$xmin)

# Need to do this bit of code before adding the next heatmap
# See - See "7.3.1 Visualize tree with multiple associated matrix" https://yulab-smu.top/treedata-book/chapter7.html
lin_hm <- lin_hm + ggnewscale::new_scale_fill() 

# Add DR status
dr_status_hm <- gheatmap(lin_hm, dr_status_data,
                         width = width,
                         offset = os,
                         color = NULL,
                         colnames_position = "top",
                         colnames_angle = angle, 
                         colnames_offset_y = 1,
                         hjust = 0,
                         font.size = font_sz) +
  scale_fill_manual(values = dr_status_colours, breaks = names(dr_status_colours) )+
  labs(fill = dr_status_lab)

# Add DST data

dr_status_hm <- dr_status_hm + ggnewscale::new_scale_fill() 

dr_data_hm <- gheatmap(dr_status_hm, dr_data, 
         # Increase offset
         offset = os*2,
         width = width, 
         low="white", high="black", color="black", 
         colnames_position = "top",
         colnames_angle = angle, 
         colnames_offset_y = 1,
         hjust = 0,
         font.size = font_sz)+
  # Define colours
  scale_fill_manual(values=c("white", "black"), 
                    labels = c("Sensitive", "Resistant", "NA"), na.value = "grey")+
  labs(fill = dr_data_lab)

# Add the mutations labels
final_tree <- dr_data_hm %<+% metadata + 
  geom_tiplab(aes(label = gene_mutation, colour = gene_mutation),
              linetype = NULL,
              align = T,
              offset = os*3.5,
              size = (n_samps*0.1)/font_sz)+
  scale_colour_discrete(guide = "none")+
  hexpand(.1, direction = 1) +
  vexpand(.1)

# Save
ggsave(outfile, final_tree)
















