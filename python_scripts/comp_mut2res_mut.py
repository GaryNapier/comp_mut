#!/usr/bin/env python

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

def get_counts(meta,samples,col):
    return dict(Counter([meta[s][col] for s in samples]))

def get_meta_proportion(meta,samples,column,targets):
    tmp = get_counts(meta,samples,column)
    target_count = sum([tmp.get(c,0) for c in targets])
    return round(target_count/sum(tmp.values()), 3)

def filter_vars(variants, mutation2sample, meta_dict, drug_of_interest):

    # Get stats for each variant
    stats_dict = defaultdict(dict)
    variants_passed = set()
    # for var in variants:
    for var in variants:

        # Get proportion or number of samples (in the case of lineage) per potential mutation
        samps = mutation2sample[var]

        if len(samps)<3: continue

        dst_proportion = get_meta_proportion(meta_dict,samps,drug_of_interest,['1'])
        sensitive_geno_proportion = get_meta_proportion(meta_dict,samps,'drtype',['Sensitive'])
        num_lins = len(set(resolve_lineages(get_counts(meta_dict,samps,'sublin'))))

        # Filter and add to dict
        if dst_proportion >= 0.5 and sensitive_geno_proportion <= 0.5 and num_lins > 1:
            variants_passed.add(var)
            stats_dict[var] = {'var': var, 
                               'gene': var[0], 
                               'mutation': var[1],
                               'n_samps' : len(samps), 
                               'dst_prop': dst_proportion, 
                               'dr_prop': sensitive_geno_proportion, 
                               'n_lins': num_lins}

    return (variants_passed, stats_dict)


def main(args):

    potential_comp_mut_file = args.potential_comp_mut_file
    metadata_file = args.metadata_file
    tbdb_file = args.tbdb_file
    drtypes_file = args.drtypes_file
    comp_mut_file = args.comp_mut_file
    tbprofiler_results_dir = args.tbprofiler_results_dir
    vars_exclude_file = args.vars_exclude_file
    potential_res_mut_outfile = args.potential_res_mut_outfile


    # FILES

    # potential_comp_mut_file = "metadata/ahpc_model_results.csv"
    # metadata_file = "../metadata/tb_data_18_02_2021.csv"
    # tbdb_file = "../tbdb/tbdb.csv"
    # drtypes_file = "../pipeline/db/dr_types.json"
    # comp_mut_file = '../pipeline/db/compensatory_mutations.csv'
    # tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/'
    # vars_exclude_file = 'metadata/var_exclude_katg_comp_mut.csv'
    # potential_res_mut_outfile = 'results/katg_potential_res_mut_stats.csv'

    # VARIABLES

    suffix = args.suffix
    genes_file = args.genes_file
    drug_of_interest = args.drug_of_interest
    comp_mut_genes = str(args.comp_mut_genes)

    # suffix = ".results.json"
    # genes = ('ahpC', 'katG', 'fabG1')
    # drug_of_interest = 'isoniazid'
    # comp_mut_genes = ('ahpC')


    # READ IN DATA

    # Read in ahpc GLM results file
    with open(potential_comp_mut_file, 'r') as f:
        potential_comp_mut_dict = csv_to_dict(f)

    # Convert to list of tuples
    potential_comp_mut_list = [(comp_mut_genes, var) for var in list(potential_comp_mut_dict)]

    # Read in metadata
    with open(metadata_file) as mf:
        meta_dict = csv_to_dict(mf)

    # Pull samples
    samples = list(meta_dict.keys())

    # Read in tbdb file
    with open(tbdb_file, 'r') as f:
        tbdb_dict = csv_to_dict_multi(f)

    # Read in DR types from json
    standardise_drtype = json.load(open(drtypes_file))

    # Get known compensatory mutations of interest
    compensatory_mutations = defaultdict(set)
    for row in csv.DictReader(open(comp_mut_file)):
        if row['Drug'] != drug_of_interest: continue
        compensatory_mutations[row['Drug']].add((row['Gene'],row['Mutation']))
        
    # Read in variants to exclude
    vars_exclude = get_vars_exclude(vars_exclude_file)

    # Read in genes file and store as set
    genes = set([l.strip().split()[0] for l in open(genes_file)])


    # Read in all the json data for (samples with) ahpC/katG only (>0.7 freq) and the metadata for those samples

    # Load mutation data using ('gene','change') as keys
    mutation2sample = defaultdict(set)
    sample2mutation = defaultdict(set)
    resistance_mutations = defaultdict(set)
    for s in tqdm(samples):
        file = "%s/%s%s" % (tbprofiler_results_dir, s, suffix)
        if os.path.isfile(file):
            data = json.load(open(file))
            # Skip mixed samps
            if ';' in data['sublin']: continue

            meta_dict[s]['drtype'] = data['drtype']
            meta_dict[s]['sublin'] = data['sublin']
            
            # MAKE SURE THE FOR LOOP BELOW IS INDENTED IN LINE WITH if os.path.isfile(file):
            # Otherwise adds sample s to mutation2sample etc
        
            for var in data['dr_variants'] + data['other_variants']:
                if var['gene'] not in genes: continue
                if var['freq'] < 0.7: continue
                if var['type']=='synonymous_variant': continue
                if (var['gene'], var['change']) in vars_exclude: continue
                key = (var['gene'],var['change'])
                mutation2sample[key].add(s)
                sample2mutation[s].add(key)
                if "drugs" in var:
                    for d in var["drugs"]:
                        if key in compensatory_mutations[d["drug"]]: continue
                        resistance_mutations[d["drug"]].add(key)

    # Classify potential compensatory mutations and filter 
    # GLM is only first step in identifying 'interesting' compensatory mutations
    # Need to check against tbprofiler results for each mutation 
    # e.g. if the mutation is lineage specific, then filter out

    potential_comp_mut_filtered, potential_comp_mut_stats = filter_vars(potential_comp_mut_list, mutation2sample, meta_dict, drug_of_interest)

    # Add the filtered potential compensatory mutations 
    # to the list of known compensatory mutations for the drug of interest

    compensatory_mutations[drug_of_interest].update(potential_comp_mut_filtered)

    potential_resistance_mutations = set()
    for s in tqdm(samples):
        # Get the comp, res and other variants for each sample in the full sample list
        comp_var = [var for var in sample2mutation[s] if var in compensatory_mutations[drug_of_interest]]
        res_var = [var for var in sample2mutation[s] if var in resistance_mutations[drug_of_interest]]
        other_vars = [var for var in sample2mutation[s] if var not in compensatory_mutations[drug_of_interest] and var not in resistance_mutations[drug_of_interest]]
        # If there is at least one comp variant and there are no (known) resistance variants
        if len(comp_var)>0 and len(res_var)==0:
            # If there are no 'other' variants print the sample and the comp variants
            if len(other_vars)==0:
                print("Sample, comp. variant")
                print(s,comp_var)
            # Store the 'other' vars as potential resistance variants
            for var in other_vars:
                potential_resistance_mutations.add(var)

    # Filter the potential resistance variants in the same way as filtering the potential comp. variants
    potential_res_mut_filtered, potential_res_mut_stats = filter_vars(potential_resistance_mutations, mutation2sample, meta_dict, drug_of_interest)

    # Write stats dictionary to file:
    with open(potential_res_mut_outfile, 'w') as f:
        writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(potential_res_mut_stats)))
        writer.writeheader()
        for row in potential_res_mut_stats:
            writer.writerow(potential_res_mut_stats[row])



parser = argparse.ArgumentParser(description='get novel potential resistance mutations from computational mutations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--potential-comp-mut-file', default = '', type = str, help = 'csv; first column = list of mutations e.g. c.-51G>A')
parser.add_argument('--metadata-file', default = '', type = str, help = 'csv of metadata')
parser.add_argument('--tbdb-file', default = '', type = str, help = 'csv from https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv')
parser.add_argument('--drtypes-file', default = '', type = str, help = 'json converting old WHO drug resistance types to new ones; https://github.com/GaryNapier/pipeline/blob/main/db/dr_types.json')
parser.add_argument('--comp-mut-file', default = '', type = str, help = 'csv of all known compensatory mutations; https://github.com/GaryNapier/pipeline/blob/main/db/compensatory_mutations.csv')
parser.add_argument('--tbprofiler-results-dir', default = '', type = str, help = 'directory of tbprofiler results containing one json per sample')
parser.add_argument('--vars-exclude-file', default = '', type = str, help = 'csv of gene,mutation to exclude. No header')
parser.add_argument('--potential-res-mut-outfile', default = '', type = str, help = 'name of output file of variants and their stats')
parser.add_argument('--suffix', default = '.results.json', type = str, help = 'suffix of json files in tbprofiler_results_dir')
parser.add_argument('--genes-file', default = '', type = str, help = 'file of genes of interest e.g. katG\n ahpC\n fabG1"')
parser.add_argument('--drug-of-interest', default = '', type = str, help = 'drug of interest to mutations e.g. "isoniazid"')
parser.add_argument('--comp-mut-genes', default = '', type = str, help = 'gene in which compensatory mutations lie e.g. "ahpC"')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
