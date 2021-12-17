# ahpc_mutations_regressions_data.py

# Get data for novel ahpc mutations regressions

# For each mutation, generate binary data for regression models:

# mutation 1
#         Has mutation?    INH_DST
# samp 1  0                   0
# samp 2  1                   1
# samp 3  0                   0

import json
from collections import defaultdict, Counter
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
import tbprofiler
from csv import DictReader

# def main(args):

# mutations_file = args.mutations_file
# metadata_file = args.metadata_file
# id_key = args.id_key

mutations_file = "metadata/novel_ahpc_mutations.txt"
metadata_file = "../metadata/tb_data_18_02_2021.csv"
id_key = "wgs_id"
mutaions_key = "change" # column name of mutation names e.g. "c.-101A>G"

# -------------
# READ IN DATA
# -------------

# Read in mutations data - make the mutations the key 
mutations_dict = {}
with open(mutations_file, 'r') as mf:
    mutations_reader = DictReader(mf, delimiter='\t')
    for row in mutations_reader:
        print(row)
        # Make the mutation the key, but also recapitulate the mutation in the key-values by including everything
        mutations_dict[row[mutaions_key]] = row

# Read in metadata
meta_dict = {}
with open(metadata_file) as mf:
    metadata_reader = csv.DictReader(mf)
    for row in metadata_reader:
        # Make the id the key, but also recapitulate the id in the key-values by including everything
        meta_dict[row[id_key]] = row





parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata-file', default = '', type = str, help = 'metadata file name with column of sample ids - only processes samples in this file')
parser.add_argument('--mutations-file', default = '', type = str, help = 'txt file of novel mutations found with find_ahpc_mutations.py; columns: wgs_id, drtype, lineage, gene, change, freq, inh_dst)
parser.add_argument('--id-key', default = '', type = str, help = 'column name in metadata file with sample ids')
parser.add_argument('--outfile',default="ahpc_mutations_stats.txt",type=str,help='name of output file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)








