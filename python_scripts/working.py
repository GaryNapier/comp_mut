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
import re


def flat_list(mylist):
    return [item for sublist in mylist for item in sublist]

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

def csv_to_dict(file):

    # If csv has unique rows for each entry of interest (e.g. sample), use this function

    # e.g. 

    # sample, lineage, DR_status
    # S1, lin1, XDR
    # S2, lin3, MDR
    # S3, lin4, sensitive
    # ->
    # {'S1': {'sample':'S1', 'lineage':'lin1', 'DR_status':'XDR'}, 
    # 'S2': {'sample':'S2', 'lineage':'lin3', 'DR_status':'MDR'},
    # 'S3': {'sample':'S3', 'lineage':'lin4', 'DR_status':'sensitive'}}

    reader = csv.DictReader(file)
    key = reader.fieldnames[0]
    out_dict = {}
    for row in reader:
        # Make the id the key, but also put the id in the key-values by including everything
        out_dict[row[key]] = row
    return out_dict

def csv_to_dict_multi(file):
    # If the csv file contains multiple rows for the key of interest ('Gene' in the example below), use this function
    
    # Convert from e.g.:

    # Gene,Mutation,Drug,Confers,Interaction,Literature
    # gyrB,p.Glu540Asp,moxifloxacin,resistance,,10.1128/AAC.00825-17;10.1128/JCM.06860-11
    # pncA,p.Thr100Ile,pyrazinamide,resistance,,10.1128/JCM.01214-17
    # pncA,p.Thr160Ala,pyrazinamide,resistance,,10.1128/JCM.01214-17
    # gid,p.Gly73Ala,streptomycin,resistance,,10.1128/AAC.02175-18
    # gid,p.Leu79Ser,streptomycin,resistance,,10.1128/AAC.02175-18
    # gid,p.Leu108Arg,streptomycin,resistance,,10.1128/AAC.02175-18

    # to:

    # 'rpoC': [{'Gene': 'rpoC',
    #    'Mutation': 'p.Asp485Asn',
    #    'Drug': 'rifampicin',
    #    'Confers': 'resistance',
    #    'Interaction': '',
    #    'Literature': ''},
    #   {'Gene': 'rpoC',
    #    'Mutation': 'p.Asp735Asn',
    #    'Drug': 'rifampicin',
    #    'Confers': 'resistance',
    #    'Interaction': '',
    #    'Literature': ''},...
    # 'rpsL': [{'Gene': 'rpsL',
    #    'Mutation': 'p.Arg86Pro',
    #    'Drug': 'streptomycin',
    #    'Confers': 'resistance',
    #    'Interaction': '',
    #    'Literature': ''},
    #   {'Gene': 'rpsL',
    #    'Mutation': 'p.Arg86Trp',
    #    'Drug': 'streptomycin',
    #    'Confers': 'resistance',
    #    'Interaction': '',
    #    'Literature': ''},... etc

    reader = csv.DictReader(file)
    key = reader.fieldnames[0]
    out_dict = {}
    for row in reader:
        if row[key] in out_dict.keys():
            out_dict[row[key]].append(row)
        else:
            out_dict[row[key]] = []
            out_dict[row[key]].append(row)
    return out_dict

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

def get_counts(in_dict, data_key):
    data_dict = {k:[] for k in in_dict.keys()}
    for mut in in_dict:
        for samp in in_dict[mut]:
            data_dict[mut].append(ahpc_dict[samp]['metadata'][data_key])

    data_counts = {k:[] for k in in_dict.keys()}
    for mut in data_dict:
        data_counts[mut] = dict(Counter(data_dict[mut]))
    return data_counts

def get_unique_mutations(in_dict, gene):
    # Pull all mutations from a certain gene and take unique
    mutations_set = []
    for samp in in_dict:
        mut_list = in_dict[samp]['mutations']
        mutations_set.append(mut['change'] for mut in mut_list if mut['gene'] == gene)
    # Convert to flat list 
    mutations_set = flat_list(mutations_set)
    # Sort
    mutations_set.sort()
    # Make a table (before converting to a set)
    table = Counter(mutations_set)
    # Get unique
    mutations_set = set(mutations_set)
    return (table, mutations_set)

def is_dict(in_dict):
    if not isinstance(in_dict, dict):
        raise Exception('is_dict() function: input is not a dictionary') 

def dr_filter(in_dict):

    # Test the proportion of sensitive DR types to all other DR types
    # If the majority are sensitive, then 'reject' (0) or 'accept' (1) 
    # Take a dictionary of DR type counts like this:
    # {'MDR': 15, 'Pre-MDR': 20, 'Pre-XDR': 8, 'XDR': 1, 'Other': 1}
    # If only 'Sensitive' in the dict the return accept = 0 
    # If 'Sensitive' is not in the counts then accept = 1
    # Else test the proportions

    is_dict(in_dict)

    if all('Sensitive' in x for x in in_dict):
        accept = 0
        return accept

    if 'Sensitive' not in in_dict:
        accept = 1
        return accept

    sum_non_sens = sum(in_dict[item] for item in in_dict if item != 'Sensitive')

    proportion = in_dict['Sensitive'] / (in_dict['Sensitive'] + sum_non_sens)

    if proportion >= 0.5: 
        accept = 0 
    else:
        accept = 1

    return accept

def lin_filter(in_dict):
    is_dict(in_dict)
    lin_filter_dict = {}
    for mutation in in_dict:
        if len(ahpc_lin_counts[mutation]) > 1:
            lin_filter_dict[mutation] = 1
        else:
            lin_filter_dict[mutation] = 0
    return lin_filter_dict

def dst_filter(in_dict):
    is_dict(in_dict)
    proportion = in_dict['1'] / (in_dict['0'] + in_dict['1'])

    if proportion < 0.5:
        accept = 0
    else:
        accept = 1

    return accept

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
# metadata_id_key = "wgs_id"
suffix = ".results.json"

# -------------
# READ IN DATA
# -------------

# Read in ahpc GLM results file
with open(ahpc_glm_results_file, 'r') as f:
    ahpc_glm_dict = csv_to_dict(f)


# Read in metadata
with open(metadata_file) as mf:
    meta_dict = csv_to_dict(mf)

# Read in tbdb file
with open(tbdb_file, 'r') as f:
    tbdb_dict = csv_to_dict_multi(f)

# Read in list of lineage-specific mutations (Fst = 1)
# See - https://www.codegrepper.com/code-examples/python/how+to+read+a+csv+file+from+a+url+with+python
with closing(requests.get(fst_results_url, stream=True)) as r:
    f = (line.decode('utf-8') for line in r.iter_lines())
    fst_dict = csv_to_dict_multi(f)

# Read in all the json data for (samples with) ahpC/katG only (>0.7 freq) and the metadata for those samples
all_data = {}
for samp in tqdm(meta_dict):
    tmp = []
    # Open the json file for the sample
    tmp_data = json.load(open(pp.filecheck("%s/%s%s" % (tbprofiler_results_dir, samp, suffix))))
    for var in tmp_data["dr_variants"] + tmp_data["other_variants"]:
        if var['gene'] not in ('ahpC', 'katG'): continue
        if var['freq'] < 0.7: continue
        tmp.append(var)
    # If sample meets criteria above, then tmp list will not be empty
    if len(tmp) > 0:
        # Create empty dict for next two entries
        all_data[samp] = {}
        # Append metadata dict
        all_data[samp]['metadata'] = {'wgs_id':samp,
    'inh_dst':meta_dict[samp]['isoniazid'],
    'main_lin': tmp_data['main_lin'],
    'sublin':tmp_data['sublin'],
    'country_code':meta_dict[samp]['country_code'],
    'drtype':tmp_data['drtype']}
        # Append mutations dict
        all_data[samp]['mutations'] = tmp


# -------------
# Wrangle data 
# -------------

# Find KNOWN ahpC mutations from tbdb
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
# for samp in tqdm(meta_dict):
for samp in tqdm(all_data):

    data = all_data[samp]
    tmp_ahpc_mutations = []
    for var in data['mutations']:

        # Pull ahpC: mutation is known or from GLM list (previously unknown)
        if var['gene'] != 'ahpC': continue
        if var['change'] not in all_ahpc_list: continue
        # Set drugs value if unknown
        if 'drugs' not in var: var['drugs'] = 'unknown'
        
        # Append the mutations
        tmp_ahpc_mutations.append(var)
        
    if len(tmp_ahpc_mutations) > 0:
        ahpc_dict[samp] = {}
        ahpc_dict[samp]['metadata'] = all_data[samp]['metadata']
        ahpc_dict[samp]['mutations'] = tmp_ahpc_mutations

# len - 563







# ---------------------------------------------------------------------
# Classify UNKNOWN ahpC and filter 
# GLM is only first step in identifying 'interesting' ahpC mutations
# Need to check against tbprofiler results for each mutation 
# e.g. if the mutation is lineage specific, then filter out
# ---------------------------------------------------------------------

# Make a dictionary with the unknown mutations as keys and the samples as values
unknown_ahpc_samps_dict = {k:[] for k in unknown_ahpc}
for samp in ahpc_dict:
    for var in ahpc_dict[samp]['mutations']:
        if var['change'] in unknown_ahpc:
            unknown_ahpc_samps_dict[var['change']].append(samp)

# DR types - only or mainly from sensitive 
ahpc_drtype_counts = get_counts(unknown_ahpc_samps_dict, 'drtype')

# Lineage counts - mutation is only from one (or two) lineages
ahpc_lin_counts = get_counts(unknown_ahpc_samps_dict, 'main_lin')

# DST counts
ahpc_dst_counts = get_counts(unknown_ahpc_samps_dict, 'inh_dst')







# FILTERS


# DR type filter
ahpc_dr_filter = {}
for mutation in ahpc_drtype_counts:
    ahpc_dr_filter[mutation] = dr_filter(ahpc_drtype_counts[mutation])

# Lineage filter
ahpc_lin_filter = lin_filter(ahpc_lin_counts)

# DST filter
ahpc_dst_filter = {}
for mutation in ahpc_drtype_counts:
    ahpc_dst_filter[mutation] = dst_filter(ahpc_dst_counts[mutation])




# Put together and add up scores







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







# -----------------------------------------
# Process katG mutations - get basic stats 
# -----------------------------------------

katg_table, katg_mutations_set = get_unique_mutations(ahpc_katg_dict, 'katG')


# Make a dict of samples with the katg mutations
katg_samps_dict = {k:[] for k in katg_mutations_set}
for samp in ahpc_katg_dict:
    for var in ahpc_katg_dict[samp]['mutations']:
        if var['gene'] == 'katG':
            katg_samps_dict[var['change']].append(samp)




# len = 110




# ------------------------------------------------------------------------------
# Pull all samples that have the unknown katGs (i.e. regardless of ahpC status)
# ------------------------------------------------------------------------------

all_katg = {}
for samp in tqdm(meta_dict):

    data = json.load(open(pp.filecheck("%s/%s%s" % (tbprofiler_results_dir, samp, suffix))))

    if any((var['gene'] == 'katG' and var['change'] in katg_samps_dict and var['freq'] >= 0.7) for var in data["other_variants"]):

        # Create empty list per id
        all_katg[samp] = {}

        all_katg[samp]['metadata'] = {'wgs_id':samp,
        'inh_dst':meta_dict[samp]['isoniazid'],
        'lineage':data['sublin'],
        'country_code':meta_dict[samp]['country_code'],
        'drtype':data['drtype']}

        all_katg[samp]['mutations'] = []

        for var in data["other_variants"]:

            # Pull ahpC: non-syn, >0.7 freq, mutation is known or from GLM list (previously unknown)
            if var['gene'] == 'katG' and var['change'] in katg_samps_dict and var['freq'] >= 0.7:

                # Set default value for drugs if no entry in the list of dictionaries
                var.setdefault('drugs', 'unknown')

                # Append the dictionary to the mutations list
                all_katg[samp]['mutations'].append({'gene':var['gene'],
                'change':var['change'],
                'type':var['type'],
                'freq':var['freq'], 
                'drugs':var['drugs']})
            else:
                continue
    else:
        continue



# len = 326






# parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--metadata-file', default = '', type = str, help = 'metadata file name with column of sample ids - only processes samples in this file')
# parser.add_argument('--mutations-file', default = '', type = str, help = 'txt file of novel mutations found with find_ahpc_mutations.py; columns: wgs_id, drtype, lineage, gene, change, freq, inh_dst)
# parser.add_argument('--id-key', default = '', type = str, help = 'column name in metadata file with sample ids')
# parser.add_argument('--outfile',default="ahpc_mutations_stats.txt",type=str,help='name of output file')
# parser.set_defaults(func=main)

# args = parser.parse_args()
# args.func(args)





