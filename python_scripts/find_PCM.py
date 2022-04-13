#! /usr/bin/env python 

# Find all compensatory mutations for drug of interest mutations

import json
from collections import defaultdict, Counter
import argparse
import os
import sys
import csv
import pathogenprofiler as pp
import tbprofiler
from tqdm import tqdm

# -----
# ARGS
# -----

def main(args):

    print(args)

    # tbprofiler_results_location = 'tbprofiler_pakistan_results/'
    metadata_file = args.metadata_file
    KCM_file = args.KCM_file
    drug_of_interest = args.drug_of_interest
    id_key = args.id_key
    tbprofiler_results_location = args.tbp_results
    suffix = args.suffix
    outfile = args.outfile
    # db = args.db

    # -------------
    # READ IN DATA
    # -------------

    # Read in metadata
    with open(metadata_file) as mf:
        metadata_reader = csv.DictReader(mf)
        meta_dict = {}
        for row in metadata_reader:
            # Make the id the key, but also recapitulate the id in the key-values by including everything
            meta_dict[row[id_key]] = row

    # Get known compensatory mutations data
    compensatory_mutations = defaultdict(set)
    for row in csv.DictReader(open(KCM_file)):
        if row['Drug'] != drug_of_interest: continue
        compensatory_mutations[row['Drug']].add((row['Gene'],row['Mutation']))


    # Wrangle data

    KCM = compensatory_mutations[drug_of_interest]

    # Get all genes for comp mutations for drug of interest
    genes = set([var[0] for var in KCM])

    # Pull novel comp mutations
    mutations_dict = defaultdict(dict)
    # for file in json_files_list:
    for samp in tqdm(meta_dict):
        # Open the json file for the sample
        file = "%s/%s%s" % (tbprofiler_results_location, samp, suffix)
        if os.path.isfile(file):
            data = json.load(open(file))
            dst = meta_dict[samp][drug_of_interest]
            lins = [lin['lin'] for lin in data['lineage']]
            lin = lins[len(lins) - 1] # I hate Python!
            
            for var in data["dr_variants"] + data["other_variants"]:
                # Save mutation if:
                # in genes for drug of interest, 
                # is not in KCM list,
                # is non-synonymous, 
                # does NOT already have a known drug association 
                # and is >0.7 freq
                if var["gene"] in genes \
                    and (var["gene"], var["change"]) not in KCM \
                        and var['type'] != 'synonymous_variant' \
                            and "drugs" not in var \
                                and var["freq"] >= 0.7:

                    mutations_dict[samp] = {'wgs_id': samp,
                    'drtype':data["drtype"],
                    'lineage':lin,
                    'gene':var["gene"],
                    'change':var["change"],
                    'freq':var["freq"], 
                    'drugs': drug_of_interest,   
                    'dst': dst}

    fieldnames = tuple(next(iter(mutations_dict.values())).keys())

    with open(outfile, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        # Loop over the dictionaries, appending each dictionary as a row in the file
        for id in mutations_dict:
            writer.writerow(mutations_dict[id])

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata-file', default = '', type = str, help = 'metadata file name with column of sample ids - only processes samples in this file')
parser.add_argument('--KCM-file', default = '', type = str, help = 'csv of all known compensatory mutations; https://github.com/GaryNapier/pipeline/blob/main/db/compensatory_mutations.csv')
parser.add_argument('--drug-of-interest', default = '', type = str, help = 'drug associated with the compensatory mutations')
parser.add_argument('--id-key', default = '', type = str, help = 'column name in metadata file with sample ids')
parser.add_argument('--tbp-results', default="results/",type=str,help='tbprofiler results directory (json files)')
# parser.add_argument('--db',default="tbdb",type=str,help='prefix to bed file for locus-DR associations')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.add_argument('--outfile',default="ahpc_variants.txt",type=str,help='name of output file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)


