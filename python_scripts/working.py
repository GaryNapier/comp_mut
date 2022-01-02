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
from collections import Counter

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

# Make list of katG mutationa to exclude
# "katG Arg463Leu was excluded from this calculation since it has been reported as not associated with isoniazid resistance" - https://www.nature.com/articles/s41598-018-21378-x
# -> The Susceptibility of Mycobacterium tuberculosis to Isoniazid and the Arg->Leu Mutation at Codon 463 of katG Are Not Associated - https://journals.asm.org/doi/10.1128/JCM.39.4.1591-1594.2001

katg_exclude = ["p.Arg463Leu"]

# ---------------------------------------------------------------------
# Pull unknown mutations from samples with ahpC mutations
# ---------------------------------------------------------------------

# Get samples with ahpC mutations
ahpc_dict = defaultdict(dict)
for samp in tqdm(meta_dict):

    # Open the json file for the sample
    data = json.load(open(pp.filecheck("%s/%s%s" % (tbprofiler_results_dir, samp, suffix))))

    # Including 'any()' statement to stop loop creating dictionary entry of empty lists
    # Only include samples which have at least one: 
    # ahpc, > 0.7 freq, and at least one ahpc mutation is either DR-associated (known) or one of the mutations from the GLM results
    if any((var['gene'] == 'ahpC' and var['freq'] >= 0.7 and var['change'] in all_ahpc_list) for var in data["dr_variants"] + data["other_variants"]):

        # Get DST data for sample
        inh_dst = meta_dict[samp]['isoniazid']

        # Create empty list per id
        ahpc_dict[samp] = {}

        ahpc_dict[samp]['metadata'] = {'wgs_id':samp,
        'inh_dst':inh_dst,
        'lineage':data['sublin'],
        'drtype':data['drtype']}

        ahpc_dict[samp]['mutations'] = []

        for var in data["dr_variants"] + data["other_variants"]:

            # Pull ahpC: non-syn, >0.7 freq, mutation is known or from GLM list (previously unknown)
            if var['gene'] == 'ahpC' and var['change'] in all_ahpc_list and var['freq'] >= 0.7:

                # Set default value for drugs if no entry in the list of dictionaries
                var.setdefault('drugs', 'unknown')

                # Append the dictionary to the mutations list
                ahpc_dict[samp]['mutations'].append({'gene':var['gene'],
                'change':var['change'],
                'type':var['type'],
                'freq':var['freq'], 
                'drugs':var['drugs']})
            else:
                continue
    else:
        continue

# Loop over samples with ahpC mutations and get their unknown katG mutations. Skip those with none. 
ahpc_katg_dict = defaultdict(dict)
for samp in ahpc_dict:
    # Open the json file for the sample
    data = json.load(open(pp.filecheck("%s/%s%s" % (tbprofiler_results_dir, samp, suffix))))

    # Only continue the loop if there is at least one meeting the criteria
    if any(var['gene'] == 'katG' and var['type'] != 'synonymous' and var['change'] not in katg_exclude and var['freq'] >= 0.7 and 'drugs' not in var.keys() for var in data["dr_variants"] + data["other_variants"]):

        for var in data["dr_variants"] + data["other_variants"]:
            # Pull katG mutaions: non-syn, >0.7 freq, mutation is unknown
            if var['gene'] == 'katG' and var['type'] != 'synonymous' and var['change'] not in katg_exclude and var['freq'] >= 0.7 and 'drugs' not in var.keys():
        
                # Set default value for drugs if no entry in the list of dictionaries
                var.setdefault('drugs', 'unknown')

                # Append the dictionary to the list of ahpC mutations
                ahpc_dict[samp]['mutations'].append({'gene':var['gene'],
                'change':var['change'],
                'type':var['type'],
                'freq':var['freq'], 
                'drugs':var['drugs']})
            else:
                continue
        # Add the whole dict entry to the new dict
        ahpc_katg_dict[samp] = ahpc_dict[samp]
    else:
        continue

# Data structure is now:
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
#                'drugs': 'unknown'}]}, 
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
#                'drugs': 'unknown' }]}}

# ------------------
# Process mutations
# ------------------

# Get set of unique katG mutations
katg_mutations_set = []
for samp in ahpc_katg_dict:
    mut_list = ahpc_katg_dict[samp]['mutations']
    # katg_mutations_set.append([mut['change'] for mut in mut_list if mut['gene'] == 'katG'])
    katg_mutations_set.append([mut['change'] for mut in mut_list if mut['gene'] == 'katG'])

# CHECK
# x = []
# for mut in katg_mutations_set:
#     x.append(len(mut))
# sum(x) = 205

# Pull out the sublists from the list (God-damn you Python...)
katg_mutations_set = [item for sublist in katg_mutations_set for item in sublist]
# CHECK - len(katg_mutations_set) = 205

# Sort
katg_mutations_set.sort()

# Make a table of the katG mutations
katg_table = Counter(katg_mutations_set)

# Get unique
katg_mutations_set = set(katg_mutations_set)

# len(katg_mutations_set) = 143

# Make a dict of samples with the katg mutations
katg_samps_dict = {k:[] for k in katg_mutations_set}
for samp in ahpc_katg_dict:
    for var in ahpc_katg_dict[samp]['mutations']:
        if var['gene'] == 'katG':
            katg_samps_dict[var['change']].append(samp)

        






# parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--metadata-file', default = '', type = str, help = 'metadata file name with column of sample ids - only processes samples in this file')
# parser.add_argument('--mutations-file', default = '', type = str, help = 'txt file of novel mutations found with find_ahpc_mutations.py; columns: wgs_id, drtype, lineage, gene, change, freq, inh_dst)
# parser.add_argument('--id-key', default = '', type = str, help = 'column name in metadata file with sample ids')
# parser.add_argument('--outfile',default="ahpc_mutations_stats.txt",type=str,help='name of output file')
# parser.set_defaults(func=main)

# args = parser.parse_args()
# args.func(args)





