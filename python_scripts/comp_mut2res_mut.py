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

    PCM_file = args.PCM_file
    metadata_file = args.metadata_file
    tbdb_file = args.tbdb_file
    drtypes_file = args.drtypes_file
    KCM_file = args.KCM_file
    tbprofiler_results_dir = args.tbprofiler_results_dir
    vars_exclude_file = args.vars_exclude_file
    PRM_stats_file = args.PRM_stats_file
    PRM_samples_file = args.PRM_samples_file
    binary_table_file = args.binary_table_file

    # FILES

    # potential_comp_mut_file = "results/isoniazid_PCM_results.csv"
    # metadata_file = "../metadata/tb_data_18_02_2021.csv"
    # tbdb_file = "../tbdb/tbdb.csv"
    # drtypes_file = "../pipeline/db/dr_types.json"
    # known_comp_mut_file = '../pipeline/db/compensatory_mutations.csv'
    # tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/'
    # vars_exclude_file = 'metadata/var_exclude_comp_mut.csv'
    # PRM_stats_file = 'results/PRM_stats.csv'
    # PRM_samples_file = 'results/PRM_samps.csv'
    # binary_table_file = 'results/'+drug_of_interest+'_binary_table.csv'

    # VARIABLES

    suffix = args.suffix
    drug_of_interest = args.drug_of_interest

    # READ IN DATA

    # # Read in GLM results file
    # with open(potential_comp_mut_file, 'r') as f:
    #     potential_comp_mut_dict = csv_to_dict(f)

    # Read in potential CM results file
    with open(PCM_file, 'r') as f:
        reader = csv.DictReader(f)        
        PCM_list = []
        for row in reader:
            PCM_list.append((row['gene'], row['change']))

    # Convert to list of tuples
    # potential_comp_mut_list = [(potential_comp_mut_dict[var]['gene'], potential_comp_mut_dict[var]['change']) for var in list(potential_comp_mut_dict)]

    # Read in metadata
    with open(metadata_file) as mf:
        meta_dict = csv_to_dict(mf)
        
    # Pull samples
    samples = list(meta_dict.keys())

    # # Read in tbdb file
    with open(tbdb_file, 'r') as f:
        tbdb_dict = csv_to_dict_multi(f, 'Drug')

    # Known resistance mutations
    # krm = [(var['Gene'], var['Mutation']) for var in tbdb_dict[drug_of_interest]]

    # Read in DR types from json
    standardise_drtype = json.load(open(drtypes_file))

    # Get known compensatory mutations of interest
    KCM = defaultdict(set)
    for row in csv.DictReader(open(KCM_file)):
        if row['Drug'] != drug_of_interest: continue
        KCM[row['Drug']].add((row['Gene'],row['Mutation']))
    ################################################################
    n_KCM = len(KCM[drug_of_interest])
    ################################################################

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
    genes = genes.union(set([var[0] for var in KCM[drug_of_interest]]))
    genes = genes.union(set([var[0] for var in PCM_list]))

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
                        if key in KCM[d["drug"]] or key in PCM_list: continue
                        resistance_mutations[d["drug"]].add(key)

    # Classify potential compensatory mutations and filter 
    # GLM models (filter_novel_comp_mut.R) is only first step in identifying 'interesting' compensatory mutations
    # Need to check against tbprofiler results for each mutation 
    # e.g. if the mutation is lineage specific, then filter out
    # First remove PCMs that are in the KCM list
    PCM_list = [var for var in PCM_list if var not in KCM[drug_of_interest]]
    ########################################################
    n_PCM_before_filtering = len(PCM_list) 
    ########################################################
    PCM_filtered, PCM_stats = filter_vars(PCM_list, mutation2sample, meta_dict, drug_of_interest)

    print(" --- list of potential comp mutations after filtering --- ")
    print(PCM_filtered)
    print(" --- stats: --- ")
    print(PCM_stats)
    print()

    # Add the filtered potential compensatory mutations 
    # to the list of known compensatory mutations for the drug of interest
    # KCM[drug_of_interest] + PCM_filtered
    CM = KCM[drug_of_interest].union(PCM_filtered)
    ############################################################################
    n_CM_after_filtering = len(CM)
    ############################################################################

    # ** MAIN BIT **
    # ** Go over all the samples and get the potential resistance mutations from the presence of comp mutations ** 
    PRM = set()
    # set up sample count vectors
    samps_CM = []
    for s in tqdm(samples):
    #     # Get the comp, res and other variants for each sample in the full sample list
        comp_var = [var for var in sample2mutation[s] if var in CM]
        res_var = [var for var in sample2mutation[s] if var in resistance_mutations[drug_of_interest]]
        other_vars = [var for var in sample2mutation[s] if var not in CM\
                    and var not in resistance_mutations[drug_of_interest]]

        # If there is at least one comp variant and there are no (known) resistance variants
        if len(comp_var)>0 and len(res_var)==0:
            # If there are no 'other' variants print the sample and the comp variants
            # if len(other_vars)==0:
                # print("Sample with at least one comp. mut. but no res. mutations:")
                # print("samp:", s, "comp. mut.:", comp_var)
                
            # Store the 'other' vars as potential resistance variants
            for var in other_vars:
                PRM.add(var)

    # Filter the potential resistance variants in the same way as filtering the potential comp. variants
    PRM_filtered, PRM_stats = filter_vars(PRM, mutation2sample, meta_dict, drug_of_interest)

    # Check potential resistance mutations have been found, quit if not
    if len(PRM_filtered) == 0:
        print()
        print("NO POTENTIAL RESISTANCE MUTATIONS FOUND FOR %s; QUITTING SCRIPT" % drug_of_interest)
        print()
        sys.exit()
    else:
        print()
        print("FOUND %s POTENTIAL RESISTANCE MUTATIONS FOR %s" % (len(PRM_filtered), drug_of_interest))
        print()

    # WRITE FILES

    # Make a dict of samps and metadata for samps with the potential res. mutations

    PRM_dict = {}

    for var in PRM_filtered:
        for samp in mutation2sample[var]:
            PRM_dict[samp] = {'wgs_id': samp,
                            'drug': drug_of_interest, 
                            'gene': var[0], 
                            'mutation': var[1],
                            'gene_mutation': var[0] + '-' + var[1], 
                            'main_lineage': meta_dict[samp]['main_lineage'], 
                            'sublin':meta_dict[samp]['sublin'], 
                            'country_code': meta_dict[samp]['country_code'], 
                            'drtype': meta_dict[samp]['drtype'],
                            'dst': meta_dict[samp][drug_of_interest]}

    with open(PRM_samples_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(PRM_dict)))
        writer.writeheader()
        for samp in PRM_dict:
            writer.writerow(PRM_dict[samp])

    # Remove stupid "^M" ('carriage return') character that python stupidly puts on to the end of each line, completely messing up my pipeline!!!
    # https://unix.stackexchange.com/questions/32001/what-is-m-and-how-do-i-get-rid-of-it
    cmd = "sed -i -e 's/\r//g' %s" % PRM_samples_file
    pp.run_cmd(cmd)

    # # Do sort uniq in case has been run on same drug before 
    # cmd = f'(head -n 1 {PRM_samples_file} && tail -n +2 {PRM_samples_file} | sort | uniq) > tmp.csv && mv tmp.csv {PRM_samples_file}'
    # pp.run_cmd(cmd)

    with open(PRM_stats_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(PRM_stats)))
        writer.writeheader()
        for row in PRM_stats:
            writer.writerow(PRM_stats[row])

    # Remove stupid "^M" ('carriage return') character that python stupidly puts on to the end of each line, completely messing up my pipeline!!!
    # https://unix.stackexchange.com/questions/32001/what-is-m-and-how-do-i-get-rid-of-it
    cmd = "sed -i -e 's/\r//g' %s" % PRM_stats_file
    pp.run_cmd(cmd)

    # # Do sort uniq in case has been run on same drug before 
    # cmd = f'(head -n 1 {PRM_stats_file} && tail -n +2 {PRM_stats_file} | sort | uniq) > tmp.csv && mv tmp.csv {PRM_stats_file}'
    # pp.run_cmd(cmd)


    # ----------- BINARY TABLE SUMMARY --------------

    binary_table = {}
    # CM
    for samp in sample2mutation:
        CM_vars = [var for var in sample2mutation[samp] if var in CM]
        KRM_vars = [var for var in sample2mutation[samp] if var in resistance_mutations[drug_of_interest]]
        PRM_vars = [var for var in sample2mutation[samp] if var in PRM_filtered]
    #     if len(CM_vars) > 0:
        binary_table[samp] = {'CM': CM_vars, 'KRM': KRM_vars, 'PRM': PRM_vars}

    # vars
    # CM
    CM_in_samps = set(flat_list([binary_table[samp]['CM'] for samp in binary_table]))
    KCM_in_samps = {var for var in CM_in_samps if var in KCM[drug_of_interest]}
    PCM_in_samps = {var for var in CM_in_samps if var in PCM_filtered}
    # RM
    KRM_in_samps = set(flat_list([binary_table[samp]['KRM'] for samp in binary_table]))
    PRM_in_samps = set(flat_list([binary_table[samp]['PRM'] for samp in binary_table]))

    # samps
    samps_CM = []
    samps_CM_and_KRM = []
    samps_CM_and_no_KRM = []
    samps_CM_and_no_KRM_and_PRM = []
    samps_CM_and_no_KRM_and_no_PRM = []
    for samp in binary_table:
        if len(binary_table[samp]['CM']) > 0:
            samps_CM.append(samp)
            if len(binary_table[samp]['KRM']) > 0:
                samps_CM_and_KRM.append(samp)
            if len(binary_table[samp]['KRM']) == 0:
                samps_CM_and_no_KRM.append(samp)
                if len(binary_table[samp]['PRM']) > 0:
                    samps_CM_and_no_KRM_and_PRM.append(samp)
                if len(binary_table[samp]['PRM']) == 0:
                    samps_CM_and_no_KRM_and_no_PRM.append(samp)

    summary_dict = {'n_CM_in_list_before_filtering': n_KCM + n_PCM_before_filtering,
                    'n_KCM_in_list': n_KCM, 
                    'n_PCM_in_list_before_filtering': n_PCM_before_filtering,
                                        
                    'n_CM_in_list_after_filtering': n_CM_after_filtering,
                    'n_PCM_in_list_after_filtering': len(PCM_filtered),
                    
                    'n_CM_in_samps': len(CM_in_samps), 
                    'n_KCM_in_samps': len(KCM_in_samps), 
                    'n_PCM_in_samps': len(PCM_in_samps), 
                    
                    'n_KRM_in_samps': len(KRM_in_samps), 
                    'n_PRM_in_samps': len(PRM_in_samps),
                    
                    "n_samps_PRM": len(PRM_dict), 
                    
                    "n_samps_CM": len(samps_CM), 
                    "n_samps_CM_and_KRM": len(samps_CM_and_KRM),
                    "n_samps_CM_and_no_KRM": len(samps_CM_and_no_KRM), 
                    "n_samps_CM_and_no_KRM_and_PRM": len(samps_CM_and_no_KRM_and_PRM), 
                    "n_samps_CM_and_no_KRM_and_no_PRM": len(samps_CM_and_no_KRM_and_no_PRM)}

    print()
    print("summary for", drug_of_interest)
    for x in summary_dict:
        print(x, summary_dict[x])
    print()

    with open(binary_table_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(binary_table)))
        writer.writeheader()
        for row in binary_table:
            writer.writerow(binary_table[row])

    # ------------ RARE muations in the DR genes for the drug of interest ----------------

    # i.e. example of CM samps with no PRM.
    # However these samples will likely have a mutation in the relevant DR genes because they have CM
    # But filtered because there is only one or two, or lineage specific
    # use samps_CM_and_no_KRM_and_no_PRM

    # Pull the DR mutations from these samps
    rare_vars = []
    n_samps_rare_vars = 0
    for samp in samps_CM_and_no_KRM_and_no_PRM:
    #     rare_var_samp = [var for var in sample2mutation[samp] if var[0] in genes]
        rare_var_samp = [var for var in sample2mutation[samp] if var[0] == 'katG']
        if len(rare_var_samp) > 0:
            rare_vars.append(rare_var_samp)
            n_samps_rare_vars += 1

    # How many samps have rare katg?
    n_samps_rare_vars
    # # How many distinct rare vars in this list?
    rare_var_cnt = Counter(flat_list(rare_vars))
    n_rare_vars = len(rare_var_cnt)
    # # Table of counts
    n_n_rare_vars = Counter(rare_var_cnt.values())

    print()
    print(" --- RARE MUTATIONS IN", drug_of_interest, "---")
    print("n samps with a rare", drug_of_interest,":", n_samps_rare_vars)
    print("n distinct rare", drug_of_interest, "vars:", n_rare_vars)
    print("table of rare", drug_of_interest, "var occurrences:")
    print(n_n_rare_vars)
    print()

    # ---------- KATG POST-HOC ANALYSIS ----------

    if drug_of_interest == 'isoniazid':
        co_gene = 'fabG1'
        fabg_dr_mutations = {(var['Gene'], var['Mutation']) \
                            for var in tbdb_dict[drug_of_interest] \
                            if var['Gene'] == co_gene}

        # Co-occurrence with fabG1 in samples with potential INH res. mutations
        co_gene_variants = []
        for samp in PRM_dict:
            variants = sample2mutation[samp]
            co_gene_variants.append([v for v in variants if v in fabg_dr_mutations])

        co_gene_variants = flat_list(co_gene_variants)

        fabg_count = dict(Counter(co_gene_variants))

        print("n samples having potential new resistance mutations which also have fabG1 DR mutations:")
        print(fabg_count)
        print("proportions:")
        print({count: round(fabg_count[count]/len(PRM_dict), 3) for count in fabg_count})

        # Compare to p.Ser315Thr - get proportion of samples with p.Ser315Thr that don't have a fabG1
        katg_ser315thr = ('katG', 'p.Ser315Thr')
        ser315thr_samps = []
        Ser315Thr_fabg_samps = []
        Ser315Thr_no_fabg_samps = []
        ser315thr_comp_mut_samps = []
        for samp in sample2mutation:
            mutations = sample2mutation[samp]

            # Total n samps with Ser315Thr
            if (katg_ser315thr in mutations):
                ser315thr_samps.append(samp)

            if (katg_ser315thr in mutations) and (any(var in mutations for var in fabg_dr_mutations)):
                Ser315Thr_fabg_samps.append(samp)

            if (katg_ser315thr in mutations) and not (any(var in mutations for var in fabg_dr_mutations)):
                Ser315Thr_no_fabg_samps.append(samp)
                
            # How many of the S315 samples have comp mutations?
            if (katg_ser315thr in mutations) and (any(var in mutations for var in CM)):
                ser315thr_comp_mut_samps.append(samp)

        Ser315Thr_fabg_prop = round(len(Ser315Thr_fabg_samps) / len(ser315thr_samps), 3)

        print("total n samps with a Ser315Thr mutation: ", len(ser315thr_samps))
        print("n samples with p.Ser315Thr and a fabG1 DR mutation: ", len(Ser315Thr_fabg_samps))
        print("n samples with p.Ser315Thr and no fabG1 DR mutations: ", len(Ser315Thr_no_fabg_samps))
        print("proportion: ", Ser315Thr_fabg_prop)
        print("n samps with a p.Ser315Thr mutation and a comp mutation: ", len(ser315thr_comp_mut_samps) )
        print("proportion: ", round(len(ser315thr_comp_mut_samps) / len(ser315thr_samps), 3))

        # 2 measures - fitness cost and levels of res
        # fabg1 gives extra resistance
        # therefore how 'resistant' are the new ones?


parser = argparse.ArgumentParser(description='get novel potential resistance mutations from compensatory mutations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--drug-of-interest', default = '', type = str, help = 'drug of interest to mutations e.g. "isoniazid"')
parser.add_argument('--PCM-file', default = '', type = str, help = 'csv; first column = list of mutations e.g. c.-51G>A')
parser.add_argument('--metadata-file', default = '', type = str, help = 'csv of metadata')
parser.add_argument('--tbdb-file', default = '', type = str, help = 'csv from https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv')
parser.add_argument('--drtypes-file', default = '', type = str, help = 'json converting old WHO drug resistance types to new ones; https://github.com/GaryNapier/pipeline/blob/main/db/dr_types.json')
parser.add_argument('--KCM-file', default = '', type = str, help = 'csv of all known compensatory mutations; https://github.com/GaryNapier/pipeline/blob/main/db/compensatory_mutations.csv')
parser.add_argument('--tbprofiler-results-dir', default = '', type = str, help = 'directory of tbprofiler results containing one json per sample')
parser.add_argument('--vars-exclude-file', default = '', type = str, help = 'csv of gene,mutation to exclude. No header')
parser.add_argument('--PRM-stats-file', default = '', type = str, help = 'name of output file of variants and their stats')
parser.add_argument('--PRM-samples-file', default = '', type = str, help = 'name of output file of variants and their stats')
parser.add_argument('--binary-table-file', default = '', type = str, help = 'name of output file of samples and all their types off mutation - CM, KRM, PRM, other')
parser.add_argument('--suffix', default = '.results.json', type = str, help = 'suffix of json files in tbprofiler_results_dir')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
