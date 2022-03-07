#!/usr/bin/env python

import sys
import pathogenprofiler as pp
import json
from collections import defaultdict, Counter
import argparse
import os
from tqdm import tqdm
import csv
from collections import Counter
import requests
from contextlib import closing
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
#     variants_passed = set()
    variants_passed = defaultdict(dict)
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
                
            variants_passed[var] = {'drug': drug_of_interest, 
                                   'gene': var[0], 
                                   'mutation': var[1],
                                   'samples': samps}
            
            stats_dict[var] = {'drug': drug_of_interest, 
                               'var': var, 
                               'gene': var[0], 
                               'mutation': var[1],
                               'n_samps' : len(samps), 
                               'dst_prop': dst_proportion, 
                               'dr_prop': sensitive_geno_proportion, 
                               'n_lins': num_lins}

    return (variants_passed, stats_dict)


def main(args):

    # FILES

    potential_comp_mut_file = args.potential_comp_mut_file
    metadata_file = args.metadata_file
    tbdb_file = args.tbdb_file
    drtypes_file = args.drtypes_file
    known_comp_mut_file = args.known_comp_mut_file
    tbprofiler_results_dir = args.tbprofiler_results_dir
    vars_exclude_file = args.vars_exclude_file
    potential_res_mut_stats_file = args.potential_res_mut_stats_file
    potential_res_mut_samples_file = args.potential_res_mut_samples_file

    # FILES

    # potential_comp_mut_file = "results/isoniazid_novel_comp_mut_model_results.csv"
    # metadata_file = "../metadata/tb_data_18_02_2021.csv"
    # tbdb_file = "../tbdb/tbdb.csv"
    # drtypes_file = "../pipeline/db/dr_types.json"
    # known_comp_mut_file = '../pipeline/db/compensatory_mutations.csv'
    # tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/'
    # vars_exclude_file = 'metadata/var_exclude_comp_mut.csv'
    # potential_res_mut_stats_file = 'results/potential_res_mut_stats.csv'
    # potential_res_mut_samples_file = 'results/potential_res_mut_samps.csv'

    # VARIABLES

    suffix = args.suffix
    drug_of_interest = args.drug_of_interest

    # READ IN DATA

    # Read in GLM results file
    with open(potential_comp_mut_file, 'r') as f:
        potential_comp_mut_dict = csv_to_dict(f)

    # Convert to list of tuples
    potential_comp_mut_list = [(potential_comp_mut_dict[var]['gene'], potential_comp_mut_dict[var]['term']) for var in list(potential_comp_mut_dict)]

    # Read in metadata
    with open(metadata_file) as mf:
        meta_dict = csv_to_dict(mf)
        
    # Pull samples
    samples = list(meta_dict.keys())

    # Read in tbdb file
    with open(tbdb_file, 'r') as f:
        tbdb_dict = csv_to_dict_multi(f, 'Drug')

    # Read in DR types from json
    standardise_drtype = json.load(open(drtypes_file))

    # Get known compensatory mutations of interest
    compensatory_mutations = defaultdict(set)
    for row in csv.DictReader(open(known_comp_mut_file)):
        if row['Drug'] != drug_of_interest: continue
        compensatory_mutations[row['Drug']].add((row['Gene'],row['Mutation']))

    # # Pull unique comp mut genes
    # known_comp_mut_genes = set([var[0] for var in compensatory_mutations[drug_of_interest]])
    
    # Read in variants to exclude
    vars_exclude = get_vars_exclude(vars_exclude_file)

    # WRANGLE DATA

    # Find genes associated with drug of interest
    genes = set()
    for var in tbdb_dict[drug_of_interest]:
        genes.add(var['Gene'])

    # Concat with genes from known and potential resistance mutations
    genes = genes.union(set([var[0] for var in compensatory_mutations[drug_of_interest]]))
    genes = genes.union(set([var[0] for var in potential_comp_mut_list]))

    # Read in all the json data for (samples with) mutations in drug-of-interest genes only (>0.7 freq) and the metadata for those samples

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
            
            # Update metadata
            meta_dict[s]['drtype'] = standardise_drtype[data['drtype']]
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
                        if d["drug"] not in drug_of_interest: continue
                        if key in compensatory_mutations[d["drug"]] or key in potential_comp_mut_list: continue
                        resistance_mutations[d["drug"]].add(key)

    # Classify potential compensatory mutations and filter
    # GLM models (filter_novel_comp_mut.R) is only first step in identifying 'interesting' compensatory mutations
    # Need to check against tbprofiler results for each mutation 
    # e.g. if the mutation is lineage specific, then filter out
    potential_comp_mut_filtered, potential_comp_mut_stats = filter_vars(potential_comp_mut_list, mutation2sample, meta_dict, drug_of_interest)

    # Add the filtered potential compensatory mutations 
    # to the list of known compensatory mutations for the drug of interest
    compensatory_mutations[drug_of_interest].update(potential_comp_mut_filtered)

    # ** MAIN BIT **
    # ** Go over all the samples and get the potential resistance mutations from the presence of comp mutations ** 
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

    # Check potential resistance mutations have been found, quit if not
    print()
    print("FOUND %s POTENTIAL RESISTANCE MUTATIONS FOR %s" % (len(potential_res_mut_filtered), drug_of_interest))
    print()

    if len(potential_res_mut_filtered) == 0:
        print()
        print("NO POTENTIAL RESISTANCE MUTATIONS FOUND FOR %s; QUITTING SCRIPT" % drug_of_interest)
        print()
        sys.exit()


    # WRITE FILES

    # Write dictionary to file:
    if os.path.isfile(potential_res_mut_samples_file):
        with open(potential_res_mut_samples_file, 'a') as f:
            writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(potential_res_mut_filtered)))
            for var in potential_res_mut_filtered:
                samps = potential_res_mut_filtered[var]['samples']
                for samp in samps:
                    writer.writerow({'drug': potential_res_mut_filtered[var]['drug'],
                        'gene': potential_res_mut_filtered[var]['gene'], 
                        'mutation': potential_res_mut_filtered[var]['mutation'], 
                        'samples': samp})
            
    else:
        with open(potential_res_mut_samples_file, 'w') as f:
            writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(potential_res_mut_filtered)))
            writer.writeheader()
            for var in potential_res_mut_filtered:
                samps = potential_res_mut_filtered[var]['samples']
                for samp in samps:
                    writer.writerow({'drug': potential_res_mut_filtered[var]['drug'],
                        'gene': potential_res_mut_filtered[var]['gene'], 
                        'mutation': potential_res_mut_filtered[var]['mutation'], 
                        'samples': samp})

    # Do sort uniq in case has been run on same drug before 
    cmd = f'(head -n 1 {potential_res_mut_samples_file} && tail -n +2 {potential_res_mut_samples_file} | sort | uniq) > tmp.csv && mv tmp.csv {potential_res_mut_samples_file}'
    pp.run_cmd(cmd)

    # Write stats dictionary to file:
    if os.path.isfile(potential_res_mut_stats_file):
        with open(potential_res_mut_stats_file, 'a') as f:
            writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(potential_res_mut_stats)))
            for row in potential_res_mut_stats:
                writer.writerow(potential_res_mut_stats[row])
    else:
        with open(potential_res_mut_stats_file, 'w') as f:
            writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(potential_res_mut_stats)))
            writer.writeheader()
            for row in potential_res_mut_stats:
                writer.writerow(potential_res_mut_stats[row])

    # Do sort uniq in case has been run on same drug before 
    cmd = f'(head -n 1 {potential_res_mut_stats_file} && tail -n +2 {potential_res_mut_stats_file} | sort | uniq) > tmp.csv && mv tmp.csv {potential_res_mut_stats_file}'
    pp.run_cmd(cmd)

parser = argparse.ArgumentParser(description='get novel potential resistance mutations from compensatory mutations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--drug-of-interest', default = '', type = str, help = 'drug of interest to mutations e.g. "isoniazid"')
parser.add_argument('--potential-comp-mut-file', default = '', type = str, help = 'csv; first column = list of mutations e.g. c.-51G>A')
parser.add_argument('--metadata-file', default = '', type = str, help = 'csv of metadata')
parser.add_argument('--tbdb-file', default = '', type = str, help = 'csv from https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv')
parser.add_argument('--drtypes-file', default = '', type = str, help = 'json converting old WHO drug resistance types to new ones; https://github.com/GaryNapier/pipeline/blob/main/db/dr_types.json')
parser.add_argument('--known-comp-mut-file', default = '', type = str, help = 'csv of all known compensatory mutations; https://github.com/GaryNapier/pipeline/blob/main/db/compensatory_mutations.csv')
parser.add_argument('--tbprofiler-results-dir', default = '', type = str, help = 'directory of tbprofiler results containing one json per sample')
parser.add_argument('--vars-exclude-file', default = '', type = str, help = 'csv of gene,mutation to exclude. No header')
parser.add_argument('--potential-res-mut-stats-file', default = '', type = str, help = 'name of output file of variants and their stats')
parser.add_argument('--potential-res-mut-samples-file', default = '', type = str, help = 'name of output file of variants and their stats')
parser.add_argument('--suffix', default = '.results.json', type = str, help = 'suffix of json files in tbprofiler_results_dir')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
