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

def get_vars_exclude(vars_exclude_file, do_lineage):

    if do_lineage == 1:
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
    
    else:

        # Read in variants to be excluded
        vars_exclude = []
        for l in open(vars_exclude_file):
            vars_exclude.append(tuple(l.strip().split(',')))

    return vars_exclude

def get_counts(meta,samples,col):
    return dict(Counter([meta[s][col] for s in samples]))

def get_meta_proportion(meta,samples,column,targets):
    # Make a dict of counts for each category (e.g. {'Sensitive': 4, 'MDR-TB': 2.. etc})
    tmp = get_counts(meta,samples,column)
    # Remove NA values (especially occurring in DST)
    tmp.pop('NA', None)
    # Get the values for the target (e.g. 'Sensitive' or '1'), default to 0 if no key available
    target_count = sum([tmp.get(c,0) for c in targets])
    # Make sure no division by 0
    if len(tmp) == 0:
        return None
    else:
        # Get the proportion of target to the rest
        return round(target_count/sum(tmp.values()), 3)

def filter_vars(variants, mutation2sample, meta_dict, drug_of_interest, genes, do_lineage):

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
        if do_lineage == 1:
                
            if dst_proportion is not None:
                
                if dst_proportion >= 0.5 and sensitive_geno_proportion <= 0.5 and var[0] in genes and num_lins > 1: 
                    variants_passed.add(var)
                    stats_dict[var] = {'drug': drug_of_interest, 
                                        'gene': var[0], 
                                        'mutation': var[1],
                                        'gene_mutation': var[0] + '-' + var[1],
                                        'n_samps' : len(samps), 
                                        'dst_prop': dst_proportion, 
                                        'dr_prop': sensitive_geno_proportion, 
                                        'n_lins': num_lins}
                    
            else:
                if sensitive_geno_proportion <= 0.5 and var[0] in genes: 
                    variants_passed.add(var)
                    stats_dict[var] = {'drug': drug_of_interest, 
                                        'gene': var[0], 
                                        'mutation': var[1],
                                        'gene_mutation': var[0] + '-' + var[1],
                                        'n_samps' : len(samps), 
                                        'dst_prop': dst_proportion, 
                                        'dr_prop': sensitive_geno_proportion, 
                                        'n_lins': num_lins}

        else:
            
            if dst_proportion is not None:
                
                if dst_proportion >= 0.5 and sensitive_geno_proportion <= 0.5 and var[0] in genes:
                    variants_passed.add(var)
                    stats_dict[var] = {'drug': drug_of_interest, 
                                        'gene': var[0], 
                                        'mutation': var[1],
                                        'gene_mutation': var[0] + '-' + var[1],
                                        'n_samps' : len(samps), 
                                        'dst_prop': dst_proportion, 
                                        'dr_prop': sensitive_geno_proportion, 
                                        'n_lins': num_lins}
                        
                        
            else:
                if sensitive_geno_proportion <= 0.5 and var[0] in genes:
                    variants_passed.add(var)
                    stats_dict[var] = {'drug': drug_of_interest, 
                                        'gene': var[0], 
                                        'mutation': var[1],
                                        'gene_mutation': var[0] + '-' + var[1],
                                        'n_samps' : len(samps), 
                                        'dst_prop': dst_proportion, 
                                        'dr_prop': sensitive_geno_proportion, 
                                        'n_lins': num_lins}

    return (variants_passed, stats_dict)

def main(args):

    # FILES

    # in 
    CM_file = args.CM_file
    KCM_file = args.KCM_file
    metadata_file = args.metadata_file
    tbdb_file = args.tbdb_file
    CM_RM_assoc_file = args.CM_RM_assoc_file
    drtypes_file = args.drtypes_file
    tbprofiler_results_dir = args.tbprofiler_results_dir
    vars_exclude_file = args.vars_exclude_file
    # out
    CM_stats_file = args.CM_stats_file
    PRM_stats_file = args.PRM_stats_file
    PRM_samples_file = args.PRM_samples_file
    binary_table_file = args.binary_table_file
    # summary_file = args.summary_file

    # VARIABLES

    suffix = args.suffix
    drug_of_interest = args.drug_of_interest
    do_lineage = int(args.do_lineage)

    # READ IN DATA

    # Read in potential CM results file
    with open(CM_file, 'r') as f:
        reader = csv.DictReader(f)        
        CM_list = []
        for row in reader:
            CM_list.append((row['gene'], row['change']))

    # Read in metadata
    with open(metadata_file) as mf:
        meta_dict = csv_to_dict(mf)

    # Pull samples
    samples = list(meta_dict.keys())

    # # Read in tbdb file
    with open(tbdb_file, 'r') as f:
        tbdb_dict = csv_to_dict_multi(f, 'Drug')

    # Read in CM-RM gene association file and store genes
    CM_RM_assoc_file
    with open(CM_RM_assoc_file, 'r') as f:
        CM_RM_assoc_dict = csv_to_dict_multi(f, 'drug')

    # WRANGLE DATA

    RM_genes = set()
    CM_genes = set()
    for x in CM_RM_assoc_dict[drug_of_interest]:
        RM_genes.add(x['RM_gene'])
        CM_genes.add(x['CM_gene'])

    # Read in DR types from json
    standardise_drtype = json.load(open(drtypes_file))

    # Get known compensatory mutations of interest
    KCM = defaultdict(set)
    for row in csv.DictReader(open(KCM_file)):
        if row['Drug'] != drug_of_interest: continue
        KCM[row['Drug']].add((row['Gene'],row['Mutation']))
    
    # Read in variants to exclude
    vars_exclude = get_vars_exclude(vars_exclude_file, do_lineage)

    # Find genes associated with drug of interest
    genes = set()
    for var in tbdb_dict[drug_of_interest]:
        genes.add(var['Gene'])

    # Concat with genes from known and potential resistance mutations
    genes = genes.union(set([var[0] for var in KCM[drug_of_interest]]))
    genes = genes.union(set([var[0] for var in CM_list]))

    # Read in all the json data for (samples with) mutations in drug-of-interest genes only (>0.7 freq) and the metadata for those samples

    # Load mutation data using ('gene','change') as keys
    mutation2sample = defaultdict(set)
    sample2mutation = defaultdict(set)
    resistance_mutations = defaultdict(set)
    samples_with_files = set()
    for s in tqdm(samples):
        file = "%s/%s%s" % (tbprofiler_results_dir, s, suffix)
        if os.path.isfile(file):
            data = json.load(open(file))
            
            # Keep record of which samples have files 
            samples_with_files.add(s)
            
            # Update metadata
            meta_dict[s]['drtype'] = standardise_drtype[data['drtype']]
            meta_dict[s]['sublin'] = data['sublin']
            
            # Skip mixed samps
            if ';' in data['sublin']: continue
            # Skip animal samples
            if 'La' in data['main_lin']: continue
            # Skip lineages 7, 8, 9
            if data['main_lin'] in {'lineage7', 'lineage8', 'lineage9'}: continue

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
                        if key in CM_list: continue
                        resistance_mutations[d["drug"]].add(key)

    # Clean up samples list - needs to agree with samples which actually have files
    samples = [samp for samp in samples if samp in samples_with_files]

    # Clean up meta_dict - needs to agree with samples which actually have files
    meta_dict = {samp: meta_dict[samp] for samp in samples}

    # Filter all CM
    CM_filtered, CM_stats = filter_vars(CM_list, mutation2sample, meta_dict, drug_of_interest, CM_genes, do_lineage)

    print()
    print("USING %s COMPENSATORY MUTATIONS FOR %s" % (len(CM_filtered), drug_of_interest))
    print()

    # WRITE FILES
    with open(CM_stats_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(CM_stats)))
        writer.writeheader()
        for row in CM_stats:
            writer.writerow(CM_stats[row])

    clean_file(CM_stats_file)

    # ** MAIN BIT **
    # ** Go over all the samples and get the potential resistance mutations from the presence of comp mutations ** 
    PRM = set()
    # set up sample count vectors
    samps_CM = []
    for s in tqdm(samples):
        # Get the comp, res and other variants for each sample in the full sample list
        comp_var = [var for var in sample2mutation[s] if var in CM_filtered]
        res_var = [var for var in sample2mutation[s] if var in resistance_mutations[drug_of_interest]]
        other_vars = [var for var in sample2mutation[s] if var not in CM_filtered\
                    and var not in resistance_mutations[drug_of_interest]]

        # Make sure there are no RM genes in res_var e.g. all the res_var can be all fabG1 and no katG
        res_var_in_RM_genes = {var[0] for var in res_var if var[0] in RM_genes}
        
        # If there is at least one comp variant and there are no (known) resistance variants
        if len(comp_var)>0 and len(res_var_in_RM_genes)==0:
                
            # Store the 'other' vars as potential resistance variants
            for var in other_vars:
                # Filter out promoters
                if 'c.-' not in var[1]:
                    PRM.add(var)

    # Filter the potential resistance variants in the same way as filtering the potential comp. variants
    PRM_filtered, PRM_stats = filter_vars(PRM, mutation2sample, meta_dict, drug_of_interest, RM_genes, do_lineage)

    # Check potential resistance mutations have been found, quit if not
    if len(PRM_filtered) > 0:

        print()
        print("FOUND %s POTENTIAL RESISTANCE MUTATIONS FOR %s" % (len(PRM_filtered), drug_of_interest))
        print()

        # WRITE PRM FILES

        with open(PRM_stats_file, 'w') as f:
            writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(PRM_stats)))
            writer.writeheader()
            for row in PRM_stats:
                writer.writerow(PRM_stats[row])

        clean_file(PRM_stats_file)

    else:
        print()
        print("NO POTENTIAL RESISTANCE MUTATIONS FOUND FOR %s" % drug_of_interest)
        print()

    # ----------- BINARY TABLE SUMMARY --------------

    binary_table = {}

    # CM
    # for samp in sample2mutation:
    for samp in samples:
        CM_vars = [var for var in sample2mutation[samp] if var in CM_filtered]
        KRM_vars = [var for var in sample2mutation[samp] if var in resistance_mutations[drug_of_interest]]
        PRM_vars = [var for var in sample2mutation[samp] if var in PRM_filtered]
        # Take out promoters in RIF
        if drug_of_interest == 'rifampicin':
            other_vars = [var for var in sample2mutation[samp] \
                        if var not in CM_filtered and \
                        var not in resistance_mutations[drug_of_interest] and \
                        var not in PRM_filtered and \
                        'c.-' not in var[1]]
        else:
            other_vars = [var for var in sample2mutation[samp] \
                        if var not in CM_filtered and \
                        var not in resistance_mutations[drug_of_interest] and \
                        var not in PRM_filtered]

        binary_table[samp] = {'wgs_id': samp, \
                            'main_lineage': meta_dict[samp]['main_lineage'], \
                            'sublin': meta_dict[samp]['main_lineage'], \
                            'country_code': meta_dict[samp]['country_code'], \
                            'drtype': meta_dict[samp]['drtype'], \
                            'dst': meta_dict[samp][drug_of_interest], \
                            'CM': CM_vars, \
                            'KRM': KRM_vars, \
                            'PRM': PRM_vars, \
                            'other_vars': other_vars}

    # Get PRM samples
    samps_PRM = []
    for samp in binary_table:
        if len(binary_table[samp]['PRM']) >0:
            samps_PRM.append(samp) 

    # Split out other vars in case more than one in there
    for samp in binary_table:
        other_vars = binary_table[samp]['other_vars']
        other_var_genes = [x[0] for x in other_vars]
        for gene in other_var_genes:
            vars_gene = [var for var in other_vars if var[0] == gene]
            binary_table[samp].update({"other_"+gene: vars_gene})
            
    # Split out KRM
    for samp in binary_table:
        KRM_vars = binary_table[samp]['KRM']
        KRM_genes = [x[0] for x in KRM_vars]
        for gene in KRM_genes:
            vars_gene = [var for var in KRM_vars if var[0] == gene]
            binary_table[samp].update({"KRM_"+gene: vars_gene})

    # Fill out the keys for the rest of the samples, getting all unique keys first
    binary_table_keys = []
    for samp in binary_table:
        binary_table_keys.append(list(binary_table[samp].keys()))
    binary_table_keys = set(flat_list(binary_table_keys))

    for samp in binary_table:
        keys_current = set(binary_table[samp].keys())
        keys_needed = [key for key in binary_table_keys if key not in keys_current]
        for key in keys_needed: 
            binary_table[samp].update({key: []})

    if len(samps_PRM) > 0:
        with open(PRM_samples_file, 'w') as f:
            for samp in samps_PRM:
                f.write("%s\n" % samp)
        clean_file(PRM_samples_file)

    with open(binary_table_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(binary_table)))
        writer.writeheader()
        for row in binary_table:
            writer.writerow(binary_table[row])
    clean_file(binary_table_file)


parser = argparse.ArgumentParser(description='get novel potential resistance mutations from compensatory mutations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# in
parser.add_argument('--do-lineage', type = str, help = 'run script analysing lineage? yes = 1, 0 = no')
parser.add_argument('--suffix', default = '.results.json', type = str, help = 'suffix of json files in tbprofiler_results_dir')
parser.add_argument('--tbprofiler-results-dir', default = '', type = str, help = 'directory of tbprofiler results containing one json per sample')
parser.add_argument('--drug-of-interest', default = '', type = str, help = 'drug of interest to mutations e.g. "isoniazid"')
parser.add_argument('--CM-file', default = '', type = str, help = 'csv; first column = list of mutations e.g. c.-51G>A')
parser.add_argument('--KCM-file', default = '', type = str, help = 'csv of all known compensatory mutations; https://github.com/GaryNapier/pipeline/blob/main/db/compensatory_mutations.csv')
parser.add_argument('--metadata-file', default = '', type = str, help = 'csv of metadata')
parser.add_argument('--tbdb-file', default = '', type = str, help = 'csv from https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv')
parser.add_argument('--CM-RM-assoc-file', default = '', type = str, help = 'csv associating CM genes to RM genes by drug, e.g. isoniazid, ahpC, katG')
parser.add_argument('--drtypes-file', default = '', type = str, help = 'json converting old WHO drug resistance types to new ones; https://github.com/GaryNapier/pipeline/blob/main/db/dr_types.json')
parser.add_argument('--vars-exclude-file', default = '', type = str, help = 'csv of gene,mutation to exclude. No header')
# out
parser.add_argument('--CM-stats-file', default = '', type = str, help = 'name of output file of CMs and their stats')
parser.add_argument('--PRM-stats-file', default = '', type = str, help = 'name of output file of PRMs and their stats')
parser.add_argument('--PRM-samples-file', default = '', type = str, help = 'name of output file of samples with PRMs')
parser.add_argument('--binary-table-file', default = '', type = str, help = 'name of output file of samples and all their types of mutation - CM, KRM, PRM, other')
# parser.add_argument('--summary-file', default = '', type = str, help = 'name of output file of CM, KRM etc counts')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
