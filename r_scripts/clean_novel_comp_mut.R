#!/usr/bin/env Rscript

# clean_novel_comp_mut.R


# Description
# Merge results from TC (taane_code_2.R) and GN (filter_novel_comp_mut.R) to produce a simple list of genes and their mutations
# These are the new potential compensatory mutations to be passed to comp_mut2res_mut.py

# Arguments to script:
# tc_file - results from TC taane_code_2.R - prefix is drug of interest e.g. isoniazid_tc.txt
# gn_results_file - results from GN filter_novel_comp_mut.R prefix is drug of interest e.g. isoniazid_novel_comp_mut_model_results.csv
# outfile - location and name of output file

# Input:

# Main steps:
# Read in data
# Filter TC results to <0.05 waldp col
# Chop down both files to just the gene and mutation cols
# rbind both results and just save unique
# Save

# Output:
# > results
# gene        change
# ahpC     c.-142G>A
# ahpC c.-47_-46insT
# ahpC      c.-48G>A
# ahpC      c.-51G>A
# ahpC      c.-52C>A
# ...etc

# RUN:
# Rscript r_scripts/clean_novel_comp_mut.R --tc_file <tc_file> --gn_results_file <gn_results_file> --outfile <outfile>

# Setup ----

library(optparse)

source("https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Functions.R")

# Arguments ----

option_list = list(
  # EXAMPLE:
  # make_option(c("-t", "--template_file_name"), type="character", default=NULL,
  #             help="input template xml file", metavar="character"),

  make_option(c("-t", "--tc_file"), type="character", default=NULL,
              help="file of results for drug of interest from TC", metavar="character"),
  make_option(c("-g", "--gn_results_file"), type="character", default=NULL,
              help="file of results for drug of interest from GN (filter_novel_comp_mut.R)", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="name of outfile", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("ARGUMENTS:")
print(opt)
print("---")
print(str(opt))

# TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING 
# setwd("~/Documents/comp_mut/")
# drug_of_interest <- 'isoniazid'
# tc_file <- "results/isoniazid_tc.txt"
# gn_results_file <- "results/isoniazid_novel_comp_mut_model_results.csv"
# outfile <- "results/isoniazid_novel_comp_mut_merged.csv"
# TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING 


tc_file <- opt$tc_file
gn_results_file <- opt$gn_results_file
outfile <- opt$outfile

tc_results <- read.table(tc_file, header = T)
gn_results <- read.csv(gn_results_file, header = T)

# Wrangle
tc_results <- round_cols(tc_results)
tc_results <- tc_results[tc_results["waldp"] < 0.05, ]

tc_results <- tc_results[, c("gene2", "pos2")]
names(tc_results) <- c("gene", "change")

gn_results <- gn_results[, c("term", "gene")]
gn_results <- gn_results[, c("gene", "term")]
names(gn_results) <- c("gene", "change")

results <- unique(rbind(tc_results, gn_results))
results <- odr(results)

write.csv(results, file = outfile, quote = F, row.names = F)






















