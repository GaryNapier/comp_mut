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
from python_scripts.utils import *

# 'SAMEA2534433'


# def tb_data(filename, suffix, meta_dict, genes, drtype_dict):

#     # Read in all the json data for samples with given genes only (>0.7 freq); combine with the metadata for those samples
#     all_data = {}

#     for samp in tqdm(meta_dict):
#         tmp = []
#         file = "%s/%s%s" % (filename, samp, suffix)
#         if os.path.isfile(file):
#             # Open the json file for the sample, skip if can't find
#             tmp_data = json.load(open(file))

#         for var in tmp_data["dr_variants"] + tmp_data["other_variants"]:
#         # for var in tmp_data["dr_variants"] + tmp_data["other_variants"]:
#             if var['gene'] not in genes: continue
#             if var['freq'] < 0.7: continue
#             tmp.append(var)
            
#         # If sample meets criteria above, then tmp list will not be empty
#         if len(tmp) > 0:
#         # if len(tmp) > 0:
#             # Create empty dict for next two entries
#             all_data[samp] = {}
#             # Append metadata dict
#             all_data[samp]['metadata'] = {
#                 'wgs_id':samp,
#                 'inh_dst':meta_dict[samp]['isoniazid'],
#                 'main_lin': tmp_data['main_lin'],
#                 # 'main_lin': tmp_data['main_lin'],
#                 'sublin':tmp_data['sublin'],
#                 # 'sublin':tmp_data['sublin'],
#                 'country_code':meta_dict[samp]['country_code'],
#                 'drtype':drtype_dict[tmp_data['drtype']]
#                 # 'drtype':standardise_drtype[tmp_data['drtype']]
#                 }
#             # Append mutations dict
#             all_data[samp]['mutations'] = tmp

#     # Find mixed samples
#     mixed_samp_lin_dict = {}
#     for samp in all_data:
#         lin = all_data[samp]['metadata']['main_lin']
#         sublin = all_data[samp]['metadata']['sublin']
#         if ";" in lin+sublin:
#             mixed_samp_lin_dict[samp] = {'lin': lin, 'sublin': sublin}

#     # Remove from data
#     for samp in list(all_data):
#         if samp in mixed_samp_lin_dict:
#             del all_data[samp]

#     return all_data



def get_vars_exclude(vars_exclude_file):

    # URL below is the results of all Fst = 1 variants from https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00817-3
    fst_results_url = 'https://raw.githubusercontent.com/GaryNapier/tb-lineages/main/fst_results_clean_fst_1_for_paper.csv'
    # See https://www.codegrepper.com/code-examples/python/how+to+read+a+csv+file+from+a+url+with+python for pulling data from url
    with closing(requests.get(fst_results_url, stream=True)) as r:
        f = (line.decode('utf-8') for line in r.iter_lines())
        fst_dict = csv_to_dict_multi(f)
    
    lin_specific_variants = []
    for gene in fst_dict:
        if gene in genes:
            for var in fst_dict[gene]:
                lin_specific_variants.append( tuple( [ gene, reformat_mutations(var['aa_pos']) ] ) )

    # Read in variants to be excluded
    vars_exclude = []
    for l in open(vars_exclude_file):
        vars_exclude.append(tuple(l.strip().split(',')))

    # Concat
    vars_exclude = vars_exclude + lin_specific_variants 
    return vars_exclude

def get_counts(in_dict, ref_dict, data_key):
    data_dict = {k:[] for k in in_dict.keys()}
    for mut in in_dict:
        for samp in in_dict[mut]:
            data_dict[mut].append(ref_dict[samp][data_key])

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
        if len(in_dict[mutation]) > 1:
            lin_filter_dict[mutation] = 1
        else:
            lin_filter_dict[mutation] = 0
    return lin_filter_dict

def dst_filter(in_dict):

    if all(x == 'NA' for x in in_dict.keys()):
        accept = 1
        return accept

    if all(x == '1' for x in in_dict.keys()):
        accept = 1
        return accept

    if all(x == '0' for x in in_dict.keys()):
        accept = 0
        return accept

    # Avoid key errors - add keys if not there and only have combination of '0' & 'NA' or '1' & 'NA'
    if '0' not in in_dict.keys() and 'NA' in in_dict.keys():
        in_dict['0'] = 0

    if '1' not in in_dict.keys() and 'NA' in in_dict.keys():
        in_dict['1'] = 0

    is_dict(in_dict)
    proportion = in_dict['1'] / (in_dict['0'] + in_dict['1'])

    if proportion < 0.5:
        accept = 0
    else:
        accept = 1

    return accept

def invert_tb_dict(in_dict):
    inv_dict = {}
    for samp in in_dict:
        for var in in_dict[samp]['mutations']:
            if var['change'] not in inv_dict:
                inv_dict[var['change']] = [samp]
            else:
                inv_dict[var['change']].append(samp)
    return inv_dict

def tb_data(filename, suffix, samps_list, genes, vars_exclude):
    # Read in all the json data for samples with given genes only (>0.7 freq), non-syn
    all_data = defaultdict(list)
    meta_dict = defaultdict(dict)

    for samp in tqdm(samps_list):
        # Load file for each samp
        file = "%s/%s%s" % (filename, samp, suffix)
        if os.path.isfile(file):
            # Open the json file for the sample, skip if can't find
            tmp_data = json.load(open(file))
            # Remove mixed samps
            if ";" in tmp_data['sublin']: continue
            # Get metadata
            meta_dict[samp] = {
                'wgs_id':samp,
                # 'inh_dst':meta_dict[samp]['isoniazid'],
                'main_lin':tmp_data['main_lin'],
                'sublin':tmp_data['sublin'],
                # 'country_code':meta_dict[samp]['country_code'],
                'drtype':standardise_drtype[tmp_data['drtype']]}

            # Loop over all variants of interest
            for var in tmp_data["dr_variants"] + tmp_data["other_variants"]:
                if var['gene'] not in genes: continue
                if var['freq'] < 0.7: continue
                if var['type'] == 'synonymous_variant': continue
                # Store key as tuple of gene and mutation and append the sample. Ignore if in exclude list.
                key = (var['gene'], var['change'])
                if key in vars_exclude: continue
                all_data[key].append(samp)

    # Remove single samples
    for var in list(all_data):
        if len(all_data[var]) < 2:
            del all_data[var]

    return (all_data, meta_dict)

def merge_metadata(meta_dict, meta_dict_csv, drug_of_interest):
    # Merge metadata from metadata csv file and from json data (tb_data() function)
    # Need country code and DST for drug of interest

    for samp in meta_dict:
        meta_dict[samp]['dst'] = meta_dict_csv[samp][drug_of_interest]
        meta_dict[samp]['country_code'] = meta_dict_csv[samp]['country_code']
    return meta_dict


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
drtypes_file = "../pipeline/db/dr_types.json"
# tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/'
# tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/gatk/results'
tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/'
# metadata_id_key = "wgs_id"
suffix = ".results.json"
genes = ('ahpC', 'katG', 'fabG1')
vars_exclude_file = 'metadata/var_exclude_katg_comp_mut.csv'
drug_of_interest = 'isoniazid'

# -------------
# READ IN DATA
# -------------

# ----------------------------------------------------
# NEED TO CHANGE VAR NAMES TO SOMETHING GENERALISABLE
# ----------------------------------------------------

# Read in ahpc GLM results file
with open(ahpc_glm_results_file, 'r') as f:
    ahpc_glm_dict = csv_to_dict(f)

# Convert to list of tuples
ahpc_glm_list = [('ahpC', var) for var in list(ahpc_glm_dict)]

# ----------------------------------------------------
# NEED TO CHANGE VAR NAMES TO SOMETHING GENERALISABLE
# ----------------------------------------------------

# Read in metadata
with open(metadata_file) as mf:
    meta_dict_csv = csv_to_dict(mf)

# Pull samples
samples = list(meta_dict_csv.keys())

# Read in tbdb file
with open(tbdb_file, 'r') as f:
    tbdb_dict = csv_to_dict_multi(f)

# Read in DR types from json
standardise_drtype = json.load(open(drtypes_file))

# Read in variants to exclude
vars_exclude = get_vars_exclude(vars_exclude_file)

# Read in all the json data for (samples with) ahpC/katG only (>0.7 freq) and the metadata for those samples
all_data, meta_dict = tb_data(tbprofiler_results_dir, suffix, samples, genes, vars_exclude)



# TESTING
# file = "%s/%s%s" % (tbprofiler_results_dir, 'SAMEA2534433', suffix)
# data = json.load(open(file))

# DR TYPES - how to load?



# -------------
# Wrangle data 
# -------------

meta_dict = merge_metadata(meta_dict, meta_dict_csv, drug_of_interest)



# Find KNOWN ahpC mutations from tbdb
# known_ahpc = []
# for mut in tbdb_dict['ahpC']:
#     known_ahpc.append(mut['Mutation'])

# # Concat known ahpc and ahpc from GLM results
# unknown_ahpc = list(ahpc_glm_dict.keys())
# all_ahpc_list = unknown_ahpc + known_ahpc




# # Get list of lineage-specific katG mutations
# lin_katg = []
# for var in fst_dict['katG']:
#     lin_katg.append(reformat_mutations(var['aa_pos']))

# # Clean up
# lin_katg = [mutation for mutation in lin_katg if mutation is not None]

# # Make list of katG mutations to exclude

# # "katG Arg463Leu was excluded from this calculation since it has been reported as not associated with isoniazid resistance" - https://www.nature.com/articles/s41598-018-21378-x
# # -> The Susceptibility of Mycobacterium tuberculosis to Isoniazid and the Arg->Leu Mutation at Codon 463 of katG Are Not Associated - https://journals.asm.org/doi/10.1128/JCM.39.4.1591-1594.2001

# katg_exclude = ["p.Arg463Leu"] + lin_katg

# # MTBseq: a comprehensive pipeline for whole genome sequence analysis of Mycobacterium tuberculosis complex isolates Kohl 2018:
# # https://peerj.com/articles/5895/
# # https://github.com/ngs-fzb/MTBseq_source/blob/master/var/res/MTB_Resistance_Mediating.txt

# katg_exclude = katg_exclude + ['p.Trp300Cys']

# --------------------
# Pull ahpC mutations
# --------------------

# # Get samples with ahpC mutations
# ahpc_dict = defaultdict(dict)
# for samp in tqdm(all_data):

#     data = all_data[samp]
#     tmp_ahpc_mutations = []
#     for var in data['mutations']:

#         # Pull ahpC: mutation is known or from GLM list (previously unknown)
#         if var['gene'] != 'ahpC': continue
#         if var['change'] not in all_ahpc_list: continue
#         # Set drugs value if unknown
#         if 'drugs' not in var: var['drugs'] = 'unknown'
        
#         # Append the mutations
#         tmp_ahpc_mutations.append(var)
        
#     if len(tmp_ahpc_mutations) > 0:
#         ahpc_dict[samp] = {}
#         ahpc_dict[samp]['metadata'] = all_data[samp]['metadata']
#         ahpc_dict[samp]['mutations'] = tmp_ahpc_mutations

# len - 456




# ---------------------------------------------------------------------
# Classify UNKNOWN ahpC and filter 
# GLM is only first step in identifying 'interesting' ahpC mutations
# Need to check against tbprofiler results for each mutation 
# e.g. if the mutation is lineage specific, then filter out
# ---------------------------------------------------------------------

# Make a dictionary with the unknown mutations as keys and the samples as values
# unknown_ahpc_samps_dict = {k:[] for k in unknown_ahpc}
# for samp in ahpc_dict:
#     for var in ahpc_dict[samp]['mutations']:
#         if var['change'] in unknown_ahpc:
#             unknown_ahpc_samps_dict[var['change']].append(samp)




# Subset all mutations-samples data to just the ahpC unknown
ahpc_glm_samps = {var: all_data[var] for var in ahpc_glm_list}

# DR types - only or mainly from sensitive 
# ahpc_drtype_counts = get_counts(unknown_ahpc_samps_dict, ahpc_dict, 'drtype')
ahpc_drtype_counts = get_counts(ahpc_glm_samps, meta_dict, 'drtype')

# Lineage counts - mutation is only from one (or two) lineages
ahpc_lin_counts = get_counts(unknown_ahpc_samps_dict, ahpc_dict, 'sublin')

# Aggregate lin counts
for mut in ahpc_lin_counts:
    ahpc_lin_counts[mut] = resolve_lineages(ahpc_lin_counts[mut])

# DST counts
ahpc_dst_counts = get_counts(unknown_ahpc_samps_dict, ahpc_dict, 'inh_dst')

# FILTER UNKNOWN

# DR type filter
ahpc_dr_filter = {}
for mutation in ahpc_drtype_counts:
    ahpc_dr_filter[mutation] = dr_filter(ahpc_drtype_counts[mutation])

# Lineage filter
ahpc_lin_filter = lin_filter(ahpc_lin_counts)

# DST filter
ahpc_dst_filter = {}
for mutation in ahpc_dst_counts:
    ahpc_dst_filter[mutation] = dst_filter(ahpc_dst_counts[mutation])

# Put together and add up scores
# Remove from ahpc mutations list if 0 

filter_dict = {}
for mutation in unknown_ahpc_samps_dict:
    filter_dict[mutation] = [ahpc_dr_filter[mutation], ahpc_lin_filter[mutation], ahpc_dst_filter[mutation]]
    if sum(filter_dict[mutation]) < len(filter_dict[mutation]):
        all_ahpc_list.remove(mutation)

# Remove from dict
for samp in list(ahpc_dict):
    if all(x['change'] not in all_ahpc_list for x in ahpc_dict[samp]['mutations']):
        del ahpc_dict[samp]

# ----------------------------
# Pull unknown katG mutations
# ----------------------------

# Find samples with any known katG in the ahpC list
samps_with_known_katg = []
for samp in ahpc_dict:
    data = all_data[samp]
    if any(var['gene'] == 'katG' and var['drugs'] != 'unknown' for var in data['mutations'] if 'drugs' in var.keys()):
        samps_with_known_katg.append(samp)
        # len = 257

# Loop over samples with ahpC mutations and get their unknown katG mutations. 
ahpc_katg_dict = defaultdict(dict)
for samp in ahpc_dict:

    # Skip samples with known katG
    if samp in samps_with_known_katg:
        continue
    else:
        data = all_data[samp]

        for var in data['mutations']:

            if var['gene'] != 'katG': continue
            if var['gene'] == 'synonymous_variant': continue
            if var['change'] in katg_exclude: continue
            if var['freq'] < 0.7: continue
        
            # Set default value for drugs if no entry in the list of dictionaries
            var.setdefault('drugs', 'unknown')
            ahpc_dict[samp]['mutations'].append(var)

        # Add the whole dict entry to the new dict
        ahpc_katg_dict[samp] = ahpc_dict[samp]
        
# len(ahpc_katg_dict) = 199

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

# len = 129

# ------------------------------------------------------------------------------
# Pull all samples that have the unknown katGs (i.e. regardless of ahpC status)
# ------------------------------------------------------------------------------

all_katg = {}
for samp in tqdm(all_data):
    data = all_data[samp]
    tmp = []
    for var in data["mutations"]:

        # Pull katG: non-syn, >0.7 freq, mutation is known or from GLM list (previously unknown)
        if var['gene'] != 'katG': continue
        if var['change'] not in katg_samps_dict: continue
        if var['freq'] < 0.7: continue

        # Set default value for drugs if no entry in the list of dictionaries
        var.setdefault('drugs', 'unknown')
        tmp.append(var)
    
    if len(tmp) > 0:
        all_katg[samp] = {}
        all_katg[samp]['metadata'] = all_data[samp]['metadata']
        all_katg[samp]['mutations'] = tmp

# len = 160

# ---------------------------
# Get stats and filter katGs
# ---------------------------

all_katg_table, all_katg_mutations_set = get_unique_mutations(all_katg, 'katG')

# Make a dictionary with the unknown mutations as keys and the samples as values
# all_katg_samps_dict = {}
# for samp in all_katg:
#     for var in all_katg[samp]['mutations']:
#         if var['change'] not in all_katg_samps_dict:
#             all_katg_samps_dict[var['change']] = [samp]
#         else:
#             all_katg_samps_dict[var['change']].append(samp)

all_katg_samps_dict = invert_tb_dict(all_katg)

# DR types - only or mainly from sensitive 
katg_drtype_counts = get_counts(all_katg_samps_dict, all_katg, 'drtype')

# Lineage counts - mutation is only from one (or two) lineages
katg_lin_counts = get_counts(all_katg_samps_dict, all_katg, 'sublin')
# Aggregate lineages
for mut in katg_lin_counts:
    katg_lin_counts[mut] = resolve_lineages(katg_lin_counts[mut])

# DST counts
katg_dst_counts = get_counts(all_katg_samps_dict, all_katg, 'inh_dst')

# FILTER KATG

# DR
katg_dr_filter = {}
for mutation in katg_drtype_counts:
    katg_dr_filter[mutation] = dr_filter(katg_drtype_counts[mutation])

# Lin
katg_lin_filter = lin_filter(katg_lin_counts)

# DST
katg_dst_filter = {}
for mutation in katg_dst_counts:
    katg_dst_filter[mutation] = dst_filter(katg_dst_counts[mutation])

# Put together and add up scores
# Remove from katG mutations list if 0 

all_katg_mutations_list = list(all_katg_mutations_set)
katg_filter_dict = {}
for mutation in all_katg_samps_dict:
    katg_filter_dict[mutation] = [katg_dr_filter[mutation], katg_lin_filter[mutation], katg_dst_filter[mutation]]
    if sum(katg_filter_dict[mutation]) < len(katg_filter_dict[mutation]):
        all_katg_mutations_list.remove(mutation)

# Remove from dict if all the mutations for the sample are not in the edited list 
for samp in list(all_katg):
    if all(x['change'] not in all_katg_mutations_list for x in all_katg[samp]['mutations']):
        del all_katg[samp]

# Remove the remaining mutations that are not in the list 
# (i.e. some samples will have one mutation in the all_katg_mutations_list and one that was just removed)
for samp in all_katg:
    for var in list(all_katg[samp]['mutations']):
        if var['change'] not in all_katg_mutations_list:
            all_katg[samp]['mutations'].remove(var)

# Invert again
# all_katg_samps_dict = {}
# for samp in all_katg:
#     for var in all_katg[samp]['mutations']:
#         if var['change'] not in all_katg_samps_dict:
#             all_katg_samps_dict[var['change']] = [samp]
#         else:
#             all_katg_samps_dict[var['change']].append(samp)

all_katg_samps_dict = invert_tb_dict(all_katg)

# -------------------------
# Co-occurrence with fabG1
# -------------------------

# Attach fabG1 to katG
all_katg_fabg = all_katg
for samp in all_katg:
    for var in all_data[samp]['mutations']:
        if var['gene'] == 'fabG1':
            all_katg_fabg[samp]['mutations'].append(var)

# Get all samples with fabG1 co-occurence
fabg_samps_katg_mutations = {}
for samp in all_katg_fabg:
    if any(x['gene'] == 'fabG1' for x in all_katg_fabg[samp]['mutations']):
        fabg_samps_katg_mutations[samp] = []
        for var in all_katg_fabg[samp]['mutations']:
            if var['gene'] == 'katG':
                fabg_samps_katg_mutations[samp].append(var['change'])
# Invert
fabg_samps_katg_mutations_inv = invert_dict_listed(fabg_samps_katg_mutations)

# Which samples do NOT have a fabG1 mutation which do have those katG mutations which co-occur with a fabG1 mutation?
katg_samps_no_fabg = {}
for samp in all_katg_fabg:
    if not any(x['gene'] == 'fabG1' for x in all_katg_fabg[samp]['mutations']):
        katg_samps_no_fabg[samp] = [] 
        for var in all_katg_fabg[samp]['mutations']:
            katg_samps_no_fabg[samp].append(var['change'])
katg_samps_no_fabg_inv = invert_dict_listed(katg_samps_no_fabg)

# Counts 
fabg_samps_katg_mutations_cnt = {}
for mut in fabg_samps_katg_mutations_inv:
    fabg_samps_katg_mutations_cnt[mut] = len(fabg_samps_katg_mutations_inv[mut])

katg_samps_no_fabg_cnt = {}
for mut in katg_samps_no_fabg_inv:
    katg_samps_no_fabg_cnt[mut] = len(katg_samps_no_fabg_inv[mut])

# Proportions
katg_fabg_prop = {}
for mut in fabg_samps_katg_mutations_cnt:
    if mut in katg_samps_no_fabg_cnt:
        katg_fabg_prop[mut] = round(fabg_samps_katg_mutations_cnt[mut] / (katg_samps_no_fabg_cnt[mut] + fabg_samps_katg_mutations_cnt[mut]), 3)







# Take all_katg - aggregate samples - how many have fabG1?


for samp in all_katg:









# Compare to p.Ser315Thr - get proportion of samples with p.Ser315Thr that don't have a fabG1


Ser315Thr_fabg_samps = []
Ser315Thr_no_fabg_samps = []
for samp in all_data:
    mutations = all_data[samp]['mutations']
    if any((x['gene'] == 'katG' and x['change'] == 'p.Ser315Thr') and any(x['gene'] == 'fabG1' for x in mutations) for x in mutations):
        Ser315Thr_fabg_samps.append(samp)
    if any((x['gene'] == 'katG' and x['change'] == 'p.Ser315Thr') and not any(x['gene'] == 'fabG1' for x in mutations) for x in mutations):
        Ser315Thr_no_fabg_samps.append(samp)

Ser315Thr_fabg_prop = round(len(Ser315Thr_fabg_samps) / (len(Ser315Thr_fabg_samps) + len(Ser315Thr_no_fabg_samps)), 3)
# 0.158






# parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--metadata-file', default = '', type = str, help = 'metadata file name with column of sample ids - only processes samples in this file')
# parser.add_argument('--mutations-file', default = '', type = str, help = 'txt file of novel mutations found with find_ahpc_mutations.py; columns: wgs_id, drtype, lineage, gene, change, freq, inh_dst)
# parser.add_argument('--id-key', default = '', type = str, help = 'column name in metadata file with sample ids')
# parser.add_argument('--outfile',default="ahpc_mutations_stats.txt",type=str,help='name of output file')
# parser.set_defaults(func=main)

# args = parser.parse_args()
# args.func(args)




