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

# mutations_file = "metadata/novel_ahpc_mutations.txt"
# metadata_file = "../metadata/tb_data_18_02_2021.csv"
# id_key = "wgs_id"
# mutaions_key = "wgs_id" # column name of mutation names e.g. "c.-101A>G"

ahpc_glm_results_file = "metadata/ahpc_model_results.csv"
metadata_file = "../metadata/tb_data_18_02_2021.csv"
tbdb_file = "../tbdb/tbdb.csv"
tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/'
metadata_id_key = "wgs_id"
suffix = ".results.json"

# test_file = 'metadata/head_metadata.csv'


# -------------
# READ IN DATA
# -------------

# Read in ahpc GLM results file
ahpc_dict = {}
with open(ahpc_glm_results_file, 'r') as f:
    ahpc_reader = csv.DictReader(f)
    # Get the header name of the ahpC mutations
    ahpc_mutation_header = ahpc_reader.fieldnames[0].split(',')[0]
    for row in ahpc_reader:
        # Make the id the key, but also recapitulate the id in the key-values by including everything
        ahpc_dict[row[ahpc_mutation_header]] = row


# Read in metadata
meta_dict = {}
with open(metadata_file) as mf:
    metadata_reader = csv.DictReader(mf)
    for row in metadata_reader:
        # Make the id the key, but also recapitulate the id in the key-values by including everything
        meta_dict[row[metadata_id_key]] = row


# Read in tbdb file
tbdb_dict = {}
with open(tbdb_file, 'r') as f:
    tbdb_reader = csv.DictReader(f)
    for row in tbdb_reader:
        if row['Gene'] in tbdb_dict.keys():
            tbdb_dict[row['Gene']].append(row)
        else:
            tbdb_dict[row['Gene']] = []
            tbdb_dict[row['Gene']].append(row)


# -------------
# Wrangle data 
# -------------

# Find KNOWN ahpC mutations form tbdb
known_ahpc = []
for mut in tbdb_dict['ahpC']:
    known_ahpc.append(mut['Mutation'])

# Concat known ahpc and ahpc from GLM results
all_ahpc_list = list(ahpc_dict.keys()) + known_ahpc

# ---------------------------------------------------------------------
# Find which samples have any ahpC mutation 
# - either known from the tbdb file or novel from the ahpC GLM results
# ---------------------------------------------------------------------

ahpc_katg_dict = defaultdict(dict)
for samp in tqdm(meta_dict):

    # Open the json file for the sample
    data = json.load(open(pp.filecheck("%s/%s%s" % (tbprofiler_results_dir, samp, suffix))))

    # Test
    # data = json.load(open(pp.filecheck("%s/%s%s" % (tbprofiler_results_dir, 'SRR6046149', suffix))))
    # [x for x in data["dr_variants"] + data["other_variants"] if x['gene'] in ['ahpC','katG']]

    # Including 'any()' statement to stop loop creating dictionary entry of empty lists
    # Only include samples which have at least one: 
    #  ahpc, > 0.7 freq, and at least one ahpc mutation is either DR-associated (known) or one of the mutations from the GLM results
    if any((var['gene'] == 'ahpC' and var['freq'] >= 0.7 and var['change'] in all_ahpc_list) for var in data["dr_variants"] + data["other_variants"]):

        # Get DST data for sample
        inh_dst = meta_dict[samp]['isoniazid']

        # Create empty list per id
        ahpc_katg_dict[samp] = {}

        ahpc_katg_dict[samp]['metadata'] = {'wgs_id':samp,
        'inh_dst':inh_dst,
        'lineage':data['sublin'],
        'drtype':data['drtype']}

        ahpc_katg_dict[samp]['mutations'] = []

        for var in data["dr_variants"] + data["other_variants"]:

            # Only include:
            # ahpC: non-syn, >0.7 freq, mutation is known or from GLM list (previously unknown)
            # katG: non-syn, >0.7 freq, mutation is unknown
            # if (var['gene'] == 'ahpC' and var['change'] in all_ahpc_list and var['freq'] >= 0.7) or (var['gene'] == 'katG' and var['type'] != 'synonymous' and var['freq'] >= 0.7 and 'drugs' not in var.keys()):
            if (var['gene'] == 'ahpC' and var['change'] in all_ahpc_list and var['freq'] >= 0.7) or (var['gene'] == 'katG' and var['type'] != 'synonymous' and var['freq'] >= 0.7):

                # Set default value for drugs if no entry in the list of dictionaries
                var.setdefault('drugs', 'unknown')

                # Append the dictionary to the list
                ahpc_katg_dict[samp]['mutations'].append({'gene':var['gene'],
                'change':var['change'],
                'type':var['type'],
                'freq':var['freq'], 
                'drugs':var['drugs']})
            else:
                continue
    else:
        continue


# Data structure:
# x = {'ID_1': {'metadata': {'wgs_id':'ID_1',
#         'inh_dst':1,
#         'lineage':'lin1',
#         'drtype':'XDR'}, 
# 	'mutations': [{'gene': 'ahpC',
#                'change': 'c.-72C>T',
#                'type': 'non_coding',
#                'freq': 1.0,
#                'drugs': 'unknown'}, 
# 		{'gene': 'katG',
#                'change': 'p.Ser315Thr',
#                'type': 'missense',
#                'freq': 1.0,
#                'drugs': [{'type': 'drug',
#                  'drug': 'isoniazid',
#                  'confidence': 'high'}]}]}, 
# 'ID_2': {'metadata': {'wgs_id':'ID_1',
#         'inh_dst':1,
#         'lineage':'lin1',
#         'drtype':'XDR'}, 
# 	'mutations': [{'gene': 'ahpC',
#                'change': 'c.-72C>T',
#                'type': 'non_coding',
#                'freq': 1.0,
#                'drugs': 'unknown'}, 
# 		{'gene': 'katG',
#                'change': 'p.Ser315Thr',
#                'type': 'missense',
#                'freq': 1.0,
#                'drugs': [{'type': 'drug',
#                  'drug': 'isoniazid',
#                  'confidence': 'high'}]}]}}



# Filter to samples with unknown katG
unknown_katg_dict = {}
test = []
for samp in ahpc_katg_dict:
    if any((var['gene'] == 'katG' and var['drugs'] == 'unknown') for var in ahpc_katg_dict[samp]['mutations']):
        for var in ahpc_katg_dict[samp]['mutations']:
            if (var['gene'] == 'katG' and var['drugs'] == 'unknown') or var['gene'] == 'katG' and :














# ------------------
# Process mutations
# ------------------

# Get set of unique mutations
# mutations_set = []
# for key in mutations_dict:
#     mutations_set.append(mutations_dict[key]['change'])
# mutations_set = set(mutations_set)






# parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--metadata-file', default = '', type = str, help = 'metadata file name with column of sample ids - only processes samples in this file')
# parser.add_argument('--mutations-file', default = '', type = str, help = 'txt file of novel mutations found with find_ahpc_mutations.py; columns: wgs_id, drtype, lineage, gene, change, freq, inh_dst)
# parser.add_argument('--id-key', default = '', type = str, help = 'column name in metadata file with sample ids')
# parser.add_argument('--outfile',default="ahpc_mutations_stats.txt",type=str,help='name of output file')
# parser.set_defaults(func=main)

# args = parser.parse_args()
# args.func(args)





