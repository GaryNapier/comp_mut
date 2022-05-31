#!/usr/bin/env Rscript

# merge_PCM.R


# Description
# Merge results from TC (taane_code_2.R) and results of find_PCM.py (results/<drug_of_interest>_PCM_data.txt) to produce a simple list of genes and their mutations
# These are the new potential compensatory mutations to be passed to comp_mut2res_mut.py

# Arguments to script:
# tc_file - results from TC taane_code_2.R - prefix is drug of interest e.g. isoniazid_tc.txt
# PCM_data_file - results from find_PCM.py prefix is drug of interest e.g. isoniazid_PCM_data.txt
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
# Rscript r_scripts/clean_novel_comp_mut.R --tc_file <tc_file> --gn_results_file <PCM_data_file> --outfile <outfile>

# Setup ----

TESTING <- 0

library(optparse)

source("https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Functions.R")

if (TESTING){
  # TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING TESTING 
  setwd("~/Documents/comp_mut/")
  drug_of_interest <- 'isoniazid'
  tc_file <- paste0("metadata/", drug_of_interest, "_tc.txt")
  PCM_data_file <- paste0("results/", drug_of_interest, "_PCM_data.txt")
  KCM_file <- "../pipeline/db/compensatory_mutations.csv"
  # outfile <- "results/isoniazid_PCM_merged.csv"
}else{
  
  # Arguments ----
  
  option_list = list(
    make_option(c("-d", "--drug_of_interest"), type="character", default=NULL,
                help="drug of interest", metavar="character"),
    make_option(c("-t", "--tc_file"), type="character", default=NULL,
                help="file of results for drug of interest from TC", metavar="character"),
    make_option(c("-g", "--PCM_data_file"), type="character", default=NULL,
                help="file of results for drug of interest from find_PCM.py", metavar="character"),
    make_option(c("-c", "--KCM_file"), type="character", default=NULL,
                help="file of known compensatory mutations", metavar="character"),
    make_option(c("-o", "--outfile"), type="character", default=NULL,
                help="name of outfile", metavar="character")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  
  print("ARGUMENTS:")
  print(opt)
  print("---")
  print(str(opt))
  
  drug_of_interest <- opt$drug_of_interest
  tc_file <- opt$tc_file
  PCM_data_file <- opt$PCM_data_file
  KCM_file <- opt$KCM_file
  outfile <- opt$outfile
}


tc_results <- read.table(tc_file, header = T)
PCM_data <- read.table(PCM_data_file, header = T)
KCM <- read.csv(KCM_file, header = T)

# Wrangle
tc_results <- round_cols(tc_results)
tc_results <- tc_results[tc_results["waldp"] < 0.05, ]

tc_results <- tc_results[, c("gene2", "pos2")]
names(tc_results) <- c("gene", "change")

PCM_data <- PCM_data[, c("gene", "change")]

results <- unique(rbind(tc_results, PCM_data))
results <- odr(results)

# Get genes for drug of interest in KCM list
KCM_genes <- unique(KCM[KCM["Drug"] == drug_of_interest, "Gene"])

# Filter out any results that are not in the KCM genes list
results <- results[results[, "gene"] %in% KCM_genes, ]

# Filter ahpC mutations
if(drug_of_interest == "isoniazid"){
  
  # Remove non-promoters from ahpC PCM
  results <- results[!(grepl("p.", results$change)), ]
  
}

if (!TESTING){
  write.csv(results, file = outfile, quote = F, row.names = F)
}





















