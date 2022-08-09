#!/usr/bin/env Rscript

# comp_mut_tree.R

# Description
# Plot a newick tree file with ggtree
# ggtree tutorial - https://garynapier.github.io/ggtree/
# Install ggtree - https://bioconductor.org/packages/release/bioc/html/ggtree.html (code: BiocManager::install("ggtree"))

# Arguments to script:
# tree-file - Newick file of tree to plot
# metadata-file - Metadata file - long format e.g. 
# wgs_id     drug      gene mutation    gene_mutation    main_lineage    sublin country_code     drtype  dst
# ERR1034655 isoniazid katG p.Gly299Ser katG-p.Gly299Ser            4         4           bg     MDR-TB <NA>
# ERR1035268 isoniazid katG p.Asp419Tyr katG-p.Asp419Tyr            4     4.3.2           br     MDR-TB    1
# ERR1035337 isoniazid katG p.Tyr413Cys katG-p.Tyr413Cys            4   4.1.2.1           br     MDR-TB    1
# project-code - Project code - project code on which to subset rows of metadata e.g. 'isoniazid'
# column - column name in which project code occurs e.g. 'drug'
# outfile - path and name of saved png

# Input:
# Newick file
# Metadata file

# Main steps:
# Read in tree and metadata
# Subset metadata to the rows with the specified project code 
# Wrangle data ffor heatmap strips of lineage, DR status and DST
# Plot tree in ggtree with the heatmap strips
# Save as png

# Output:

# RUN:
# Rscript r_scripts/comp_mut_tree.R <tree-file> <metadata-file> <project-code> <column> <outfile>

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

setwd("~/Documents/comp_mut/r_scripts/")
# Setup
project_code <- "isoniazid"
# project_code_col_name <- "drug"
drug <- "isoniazid"

# Paths ---- 

results_path <- "../results/"
newick_path <- paste0(results_path, "newick/")

# Files ----

# in
tree_file <- paste0(newick_path, drug, ".treefile")
# metadata_file <- paste0("results/", project_code, "_PRM_samps.csv")
metadata_file <- paste0(results_path, drug, "_binary_table.csv")
# out
outfile <- paste0(results_path, drug, "_tree.png")

stopifnot(file.exists(tree_file))
stopifnot(file.exists(metadata_file))



# Read in data ----

tree_all_samps <- read.newick(tree_file)
metadata <- read.csv(metadata_file, header = T)
# metadata <- read.csv("results/isoniazid_binary_table.csv", header = T)

# Wrangle data ----

drug_abv_table <- data.frame(drug_full_nm = c("isoniazid", "rifampicin"),
                             drug_abv = c("INH", "RIF"))

drug_abv <- drug_abv_table[drug_abv_table["drug_full_nm"] == drug, "drug_abv"]

# Midpoint root the tree
tree_all_samps <- phytools::midpoint.root(tree_all_samps)

# Metadata
# Subset metadata to just the project code rows
# metadata <- metadata[metadata[, project_code_col_name] == project_code, ]

# Clean data ----

# Lineages - remove 'lineage' and convert to factor
metadata$main_lineage <- factor(gsub('lineage', '', metadata$main_lineage))
metadata$sublin <- factor(gsub('lineage', '', metadata$sublin))
# Convert DST to factor
metadata$dst <- as.factor(metadata$dst)

# Remove all the "[]" present in the binary table caused by python's lists
metadata <- clean_binary_table(metadata)

# Subset to just those samps with PRM
metadata <- subset(metadata, !(is.na(PRM)))

# Make ID rownames of metadata so heatmap strips work in ggtree
row.names(metadata) <- metadata$wgs_id

# Get data for separate heatmap strips
lin_data <- metadata[, "main_lineage", drop = F]
dr_status_data <- metadata[, "drtype", drop = F]
dr_data <- metadata[,'dst', drop = F]
PRM_data <- subset(metadata, PRM_bin == "present")
PRM_pivot <- reshape2::dcast(PRM_data, wgs_id ~ PRM, value.var = "PRM_bin")
row.names(PRM_pivot) <- PRM_pivot$wgs_id
PRM_pivot <- drop_cols(PRM_pivot, "wgs_id")
PRM_data <- ifelse(is.na(PRM_pivot), "0", "1")
colnames(PRM_data) <- gsub("katG-p.", "", colnames(PRM_data))

# Change col headers to match legends 
lin_data_lab <- "Lineage"
dr_status_lab <- "DR type"
dr_data_lab <- paste(drug_abv, "DST")

colnames(lin_data) <- lin_data_lab
colnames(dr_status_data) <- dr_status_lab
colnames(dr_data) <- "DST"

# Colours
alpha <- 0.9
lin_colours <- rainbow(length(unique(metadata$main_lineage)), alpha = alpha)
dr_df <- data.frame(drtype = c("Sensitive", "Pre-MDR-TB", "MDR-TB", "Pre-XDR-TB", "XDR-TB", "Other"), 
                    col = c("green1", "yellow2", "orange1", "red1", "black", "grey"))
dr_status_colours <- scales::alpha(gplots::col2hex(dr_df$col), alpha = alpha)
names(lin_colours) <- c(sort(unique(metadata$main_lineage)))
names(dr_status_colours) <- dr_df$drtype

# Trim the tree to what's in the metadata
tree_all_samps <- keep.tip(tree_all_samps, intersect(metadata$wgs_id, tree_all_samps$tip.label))
n_samps <- length(tree_all_samps$tip.label)

# Set up ggtree parameters ----

width <- 0.03
font_sz <- 3.25
line_sz <- 0.25
angle <- 90

x_lim <- c(0, 0.15)
y_lim <- c(-10, n_samps + (n_samps * 0.1))
legend_spec <- theme(legend.title = element_text(size = 9),
                     legend.text = element_text(size = 7),
                     legend.key.size = unit(0.3, "cm"))

max_dist <- castor::get_tree_span(tree_all_samps, as_edge_count=FALSE)$max_distance

# Make tree ----

ggtree_all_samps <- ggtree(tree_all_samps,
                           size = line_sz)+
  coord_cartesian(xlim = x_lim)
# scale_x_continuous(limits = x_lim)+

# Add lineage data 
lin_hm <- gheatmap(ggtree_all_samps, lin_data, 
                   color = NA,
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
                         color = NA,
                         width = width,
                         offset = os,
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
                       low="white", high="black", color=NA, 
                       colnames_position = "top",
                       colnames_angle = angle, 
                       colnames_offset_y = 1,
                       hjust = 0,
                       font.size = font_sz)+
  # Define colours
  scale_fill_manual(values=c("white", "black"), 
                    labels = c("Sensitive", "Resistant", "NA"), na.value = "grey")+
  labs(fill = "DST")

# Add binary presence/absence for PRMs
final_tree <- gheatmap(dr_data_hm, PRM_data,
# gheatmap(dr_data_hm, PRM_data, 
                       offset = os*3.5,
                       width = 0.75, 
                       color = "black",
                       low="white", high="black",
                       colnames = T,
                       colnames_position = "top",
                       colnames_angle = angle, 
                       colnames_offset_y = 1,
                       hjust = 0,
                       font.size = font_sz)+
  scale_fill_manual(values=c("white", "black"))+
  geom_treescale(y = -5)+
  hexpand(0.1)+
  vexpand(.35)+
  theme(legend.position = 'bottom')
  # annotate('text',
  #          x = 0.095, y = n_samps+6,
  #          label = "Putative resistance mutation\nabsent (white), present (black)",
  #          size = font_sz,
  #          hjust = 0)

# Save
ggsave(outfile, final_tree, width = 1100/3, height = 700/3, units = "mm")










