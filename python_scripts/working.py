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
import requests
from contextlib import closing

def reformat_mutations(x):
    aa_short2long = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln',
    'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
    'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
    'Y': 'Tyr', 'V': 'Val', '*': '*', '-': '-'
    }
    re_obj = re.search("([0-9]+)([A-Z\*])>([0-9]+)([A-Z\*])",x)
    if re_obj:
        codon_num = int(re_obj.group(1))
        ref = aa_short2long[re_obj.group(2)]
        alt = aa_short2long[re_obj.group(4)]
        return "p.%s%s%s" % (ref,codon_num,alt)
    else:
        return None

def invert_dict(d):
    inverse = dict()
    for key in d:
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse:
                # If not create a new list
                inverse[item] = [key]
            else:
                inverse[item].append(key)
    return inverse

def get_counts(in_dict, data_key):
    data_dict = {k:[] for k in in_dict.keys()}
    for mut in in_dict:
        for samp in in_dict[mut]:
            data_dict[mut].append(ahpc_dict[samp]['metadata'][data_key])

    data_counts = {k:[] for k in in_dict.keys()}
    for mut in data_dict:
        data_counts[mut] = dict(Counter(data_dict[mut]))
    return data_counts

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
fst_results_url = 'https://raw.githubusercontent.com/GaryNapier/tb-lineages/main/fst_results_clean_fst_1_for_paper.csv'
metadata_id_key = "wgs_id"
suffix = ".results.json"

# test_file = 'metadata/head_metadata.csv'


# -------------
# READ IN DATA
# -------------

# Read in ahpc GLM results file
ahpc_glm_dict = {}
with open(ahpc_glm_results_file, 'r') as f:
    ahpc_reader = csv.DictReader(f)
    # Get the header name of the ahpC mutations
    ahpc_mutation_header = ahpc_reader.fieldnames[0].split(',')[0]
    for row in ahpc_reader:
        # Make the id the key, but also recapitulate the id in the key-values by including everything
        ahpc_glm_dict[row[ahpc_mutation_header]] = row


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


# Read in list of lineage-specific mutations (Fst = 1)
# See - https://www.codegrepper.com/code-examples/python/how+to+read+a+csv+file+from+a+url+with+python
fst_dict = {}
with closing(requests.get(fst_results_url, stream=True)) as r:
    f = (line.decode('utf-8') for line in r.iter_lines())
    fst_reader = csv.DictReader(f, delimiter=',', quotechar='"')
    for row in fst_reader:
        if row['Gene'] in fst_dict.keys():
            fst_dict[row['Gene']].append(row)
        else:
            fst_dict[row['Gene']] = []
            fst_dict[row['Gene']].append(row)

# -------------
# Wrangle data 
# -------------

# Find KNOWN ahpC mutations form tbdb
known_ahpc = []
for mut in tbdb_dict['ahpC']:
    known_ahpc.append(mut['Mutation'])

# Concat known ahpc and ahpc from GLM results
unknown_ahpc = list(ahpc_glm_dict.keys())
all_ahpc_list = unknown_ahpc + known_ahpc

# Get list of lineage-specific katG mutations
lin_katg = []
for var in fst_dict['katG']:
    lin_katg.append(reformat_mutations(var['aa_pos']))

# Clean up
lin_katg = [mutation for mutation in lin_katg if mutation is not None]

# Make list of katG mutationa to exclude
# "katG Arg463Leu was excluded from this calculation since it has been reported as not associated with isoniazid resistance" - https://www.nature.com/articles/s41598-018-21378-x
# -> The Susceptibility of Mycobacterium tuberculosis to Isoniazid and the Arg->Leu Mutation at Codon 463 of katG Are Not Associated - https://journals.asm.org/doi/10.1128/JCM.39.4.1591-1594.2001

katg_exclude = ["p.Arg463Leu"] + lin_katg





# --------------------
# Pull ahpC mutations
# --------------------

# Get samples with ahpC mutations
ahpc_dict = defaultdict(dict)
for samp in tqdm(meta_dict):

    # Open the json file for the sample
    data = json.load(open(pp.filecheck("%s/%s%s" % (tbprofiler_results_dir, samp, suffix))))

    # Including 'any()' statement to stop loop creating dictionary entry of empty lists
    # Only include samples which have at least one: 
    # ahpc, > 0.7 freq, and at least one ahpc mutation is either DR-associated (known) or one of the mutations from the GLM results
    if any((var['gene'] == 'ahpC' and var['freq'] >= 0.7 and var['change'] in all_ahpc_list) for var in data["dr_variants"] + data["other_variants"]):

        # Create empty list per id
        ahpc_dict[samp] = {}

        ahpc_dict[samp]['metadata'] = {'wgs_id':samp,
        'inh_dst':meta_dict[samp]['isoniazid'],
        'lineage':data['sublin'],
        'country_code':meta_dict[samp]['country_code'],
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








# ----------------------
# Classify unknown ahpC
# ----------------------

# TB-Profiler



# Make a dictionary with the unknown mutations as keys and the samples as values
unknown_ahpc_samps_dict = {k:[] for k in unknown_ahpc}
for samp in ahpc_dict:
    for var in ahpc_dict[samp]['mutations']:
        if var['change'] in unknown_ahpc:
            unknown_ahpc_samps_dict[var['change']].append(samp)

# Invert the dict
# unknown_ahpc_samps_dict_inv = invert_dict(unknown_ahpc_samps_dict)

# DR types - only or mainly from sensitive 
ahpc_drtype_counts = get_counts(unknown_ahpc_samps_dict, 'drtype')
# Lineage counts - mutation is only from one (or two) lineages
ahpc_lin_counts = get_counts(unknown_ahpc_samps_dict, 'lineage')

# No/little katG co-occurence



















# ----------------------------
# Pull unknown katG mutations
# ----------------------------

# Find samples with any known katG in the ahpC list
samps_with_known_katg = []
for samp in ahpc_dict:
    data = json.load(open(pp.filecheck("%s/%s%s" % (tbprofiler_results_dir, samp, suffix))))
    if any(var['gene'] == 'katG' for var in data['dr_variants']):
        samps_with_known_katg.append(samp)
        # len = 301

# Loop over samples with ahpC mutations and get their unknown katG mutations. 
ahpc_katg_dict = defaultdict(dict)
for samp in ahpc_dict:

    # Skip samples with known katG
    if samp in samps_with_known_katg:
        continue
    else:
        data = json.load(open(pp.filecheck("%s/%s%s" % (tbprofiler_results_dir, samp, suffix))))

        # Only continue the loop if there is at least one meeting the katG criteria
        if any(var['gene'] == 'katG' and var['type'] != 'synonymous' and var['change'] not in katg_exclude and var['freq'] >= 0.7 for var in data['other_variants']):

            for var in data['other_variants']:
                # Pull katG mutaions: non-syn, >0.7 freq, mutation is unknown
                if var['gene'] == 'katG' and var['type'] != 'synonymous' and var['change'] not in katg_exclude and var['freq'] >= 0.7:
            
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

# len(ahpc_katg_dict) = 153


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
# 'ID_2': {'metadata': {'wgs_id':'ID_2',
#         'inh_dst':1,
#         'lineage':'lin1',
#         'drtype':'XDR'}, 
# 	'mutations': [{'gene': 'ahpC',
#                'change': 'c.-72C>T',
#                'type': 'non_coding',
#                'freq': 1.0,
#                'drugs': [isoniazid]}, 
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





