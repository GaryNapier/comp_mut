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

TESTING <- 1

library(optparse)
library(phytools)
library(ggplot2)
library(ggtree)
library(scales)
library(gplots)
library(castor)
library(ggnewscale)

source("https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Functions.R")

# if (!TESTING){

#   # Arguments ----
# 
#   option_list = list(
#     make_option(c("-t", "--tree_file"), type="character", default=NULL,
#                 help="path to tree file to be plotted", metavar="character"),
#     make_option(c("-m", "--metadata_file"), type="character", default=NULL,
#                 help="path to long-format metadata file containing columns:
#                 wgs_id (ids of samples),
#                 drug (name of drug to which the other columns of metadata correspond e.g. gene, dst),
#                 main_lineage (main lineage of each samp),
#                 drtype (Sensitive, pre-MDR, MDR etc),
#                 dst (binary col of drug susceptibility test corresponding to drug column)",
#                 metavar="character"),
#     make_option(c("-p", "--project_code"), type="character", default=NULL,
#                 help="enter project code on which to subset rows of metadata e.g. 'isoniazid'", metavar="character"),
#     make_option(c("-c", "--column"), type="character", default=NULL,
#                 help="column name in which project code occurs e.g. 'drug'", metavar="character"),
#     make_option(c("-o", "--outfile"), type="character", default=NULL,
#                 help="path and name of saved png", metavar="character")
#   );
#   
#   opt_parser = OptionParser(option_list=option_list);
#   opt = parse_args(opt_parser);
#   
#   print("ARGUMENTS:")
#   print(opt)
#   print("---")
#   print(str(opt))
#   
#   # Setup
#   project_code <- opt$project_code
#   # project_code_col_name <- opt$column
#   
#   # Files ----
#   
#   tree_file <- opt$tree_file
#   metadata_file <- opt$metadata_file
#   outfile <- opt$outfile
# 
# }else{

# TESTING TESTING TESTING TESTING TESTING TESTING 
setwd("~/Documents/comp_mut/")
# Setup
project_code <- "isoniazid"
# project_code_col_name <- "drug"
drugs <- c("isoniazid", "rifampicin")

for(drug in seq(drugs)){
  
  # Files ----
  # tree_file <- paste0("results/newick/", project_code, ".filt.val.gt.g.snps.fa.treefile")
  # # metadata_file <- paste0("results/", project_code, "_PRM_samps.csv")
  # metadata_file <- paste0("results/", project_code, "_binary_table.csv")
  # outfile <- paste0("results/", project_code, "_tree.png")
  
  tree_file <- paste0("results/newick/", drugs[drug], ".treefile")
  # metadata_file <- paste0("results/", project_code, "_PRM_samps.csv")
  metadata_file <- paste0("results/", drugs[drug], "_binary_table.csv")
  outfile <- paste0("results/", drugs[drug], "_tree.png")
  outfile_tiff <- paste0("results/", drugs[drug], "_tree.tiff")
  
  # TESTING TESTING TESTING TESTING TESTING TESTING 
  
  # }
  
  # Read in data ----
  
  tree_all_samps <- read.newick(tree_file)
  metadata <- read.csv(metadata_file, header = T)
  # metadata <- read.csv("results/isoniazid_binary_table.csv", header = T)
  
  # Wrangle data ----
  
  drug_abv_table <- data.frame(drug_full_nm = c("isoniazid", "rifampicin"),
                               drug_abv = c("INH", "RIF"))
  
  # drug_abv <- drug_abv_table[drug_abv_table["drug_full_nm"] == project_code, "drug_abv"]
  drug_abv <- drug_abv_table[drug_abv_table["drug_full_nm"] == drugs[drug], "drug_abv"]
  
  # Midpoint root the tree
  tree_all_samps <- phytools::midpoint.root(tree_all_samps)
  n_samps <- length(tree_all_samps$tip.label)
  
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
  
  # Change col headers to match legends 
  lin_data_lab <- "Lineage"
  dr_status_lab <- "DR type"
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
  
  
  if (drugs[drug] == "isoniazid"){
    
    sz <- 1.5
    
  }else if(drugs[drug] == "rifampicin"){
    
    sz <- 2
    
  }
  
  dr_data_hm
  
  
  # Add the mutations labels, make adjustments to spacing and add scale
  # For adding tip labels which are not sample IDs, see - https://yulab-smu.top/treedata-book/faq.html
  final_tree <- dr_data_hm %<+% metadata + 
    # geom_tiplab(aes(label = gene_mutation, colour = gene_mutation),
    geom_tiplab(aes(label = PRM, colour = PRM),
                linetype = NULL,
                align = T,
                # angle = 45, 
                offset = os*3.5,
                size = sz)+
    scale_colour_discrete(guide = "none")+
    # hexpand(.1, direction = 1) +
    hexpand(.3, direction = 1) +
    vexpand(.1)+
    geom_treescale(y = -5)+
    annotate('text', 
             x = 0, y = c(n_samps, n_samps-(n_samps*0.1)),
             # label = c(project_code, 
             #           sprintf("n = %s", n_samps)), 
             label = c(drugs[drug], 
                       sprintf("n = %s", n_samps)), 
             size = 7, 
             hjust = 0)+
    theme(legend.position = 'bottom')
  
  final_tree
  # Save
  ggsave(outfile, final_tree, width = 1100/3, height = 700/3, units = "mm")
  ggsave(outfile_tiff, final_tree, width = 1100/3, height = 700/3, units = "mm")
  
  # tiff(outfile_tiff, width = 1100/3, height = 700/3, units = "mm", res=300)
  # # tiff("test.tiff", units="in", width=5, height=5, res=300)
  # final_tree
  # dev.off()
  
}



tree_all_samps

getMRCA(tree_all_samps, c("SRR671799", "SRR8651577"))


id_PRM <- select(metadata, wgs_id, PRM)


id_PRM_split <- split(id_PRM, id_PRM$PRM)

lapply(id_PRM_split, function(x){
  getMRCA(tree_all_samps, x$wgs_id)
})



# $`katG-p.Met257Val; katG-p.Gln88Pro`
# NULL
# $`katG-p.Phe183Leu`
# NULL
# $`katG-p.Pro232Ser`
# [1] 178
# $`katG-p.Pro232Ser; katG-p.Glu233Gly; katG-p.Asp675Tyr`
# NULL
# $`katG-p.Thr271Ile`
# [1] 182
# $`katG-p.Thr326Pro`
# NULL







p <- ggtree(tree_all_samps, size = line_sz) + xlim(NA, 0.1)
# p <- p + geom_tiplab(align = T, size = 2)
p <- p + geom_text(aes(label=node), hjust=-.3)
p + geom_cladelab(node=190, label="test label", align = T)



# 
# # for(i in 1:2){
# #   # tiff(sprintf("test_%s.tiff", i), units="in", width=5, height=5, res=300)
# #   # ggplot(data.frame(x = 1:10, y = 1:10), aes(x = x, y = y))+geom_point()
# #   # plot(1:10)
# #   x <- ggplot(data.frame(x = 1:10, y = 1:10), aes(x = x, y = y))+geom_point()
# #   # dev.off()
# #   ggsave(sprintf("test_%s.png", i), x)
# # }
# # 
# # 
# # 
# final_tree <- dr_data_hm %<+% metadata +
#   # geom_tiplab(aes(label = gene_mutation, colour = gene_mutation),
#   geom_tiplab(aes(label = PRM, colour = PRM),
#               linetype = NULL,
#               align = T,
#               # angle = 45,
#               offset = os*3.5,
#               size = sz)+
#   scale_colour_discrete(guide = "none")+
#   # hexpand(.1, direction = 1) +
#   hexpand(.3, direction = 1) +
#   vexpand(.1)+
#   geom_treescale(y = -5)+
#   annotate('text',
#            x = 0, y = c(n_samps, n_samps-(n_samps*0.1)),
#            # label = c(project_code,
#            #           sprintf("n = %s", n_samps)),
#            label = c(drugs[drug],
#                      sprintf("n = %s", n_samps)),
#            size = 7,
#            hjust = 0)+
#   theme(legend.position = 'bottom')
# 
# 
# 
# 
# data(iris)
# rn <- paste0(iris[,5], "_", 1:150)
# rownames(iris) <- rn
# d_iris <- dist(iris[,-5], method="man")
# 
# tree_iris <- ape::bionj(d_iris)
# grp <- list(setosa     = rn[1:50],
#             versicolor = rn[51:100],
#             virginica  = rn[101:150])
# 
# p_iris <- ggtree(tree_iris)
# groupOTU(p_iris, grp, 'Species') + aes(color=Species)
# 
# 
# 
# final_tree <- dr_data_hm %<+% metadata +
#   # geom_tiplab(aes(label = gene_mutation, colour = gene_mutation),
#   geom_tiplab(aes(label = PRM, colour = PRM),
#               linetype = NULL,
#               align = T,
#               # angle = 45,
#               offset = os*3.5,
#               size = sz)+
#   scale_colour_discrete(guide = "none")+
#   # hexpand(.1, direction = 1) +
#   hexpand(.3, direction = 1) +
#   vexpand(.1)+
#   geom_treescale(y = -5)+
#   annotate('text',
#            x = 0, y = c(n_samps, n_samps-(n_samps*0.1)),
#            # label = c(project_code,
#            #           sprintf("n = %s", n_samps)),
#            label = c(drugs[drug],
#                      sprintf("n = %s", n_samps)),
#            size = 7,
#            hjust = 0)+
#   theme(legend.position = 'bottom')
# 
# 
# 
# 
# x <- dplyr::select(metadata, wgs_id, PRM)
# 
# x <- split(x, x$PRM)
# 
# x <- lapply(x, function(x){
#   x$wgs_id
# })
# 


# test <- groupOTU(final_tree, x, 'PRM') + aes(colour = PRM)






  