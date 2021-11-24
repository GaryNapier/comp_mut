#! /usr/bin/env python 

# Find all ahpC mutations

import json
from collections import defaultdict, Counter
import argparse
import os
import sys
import csv
import pathogenprofiler as pp
import tbprofiler

# -----
# ARGS
# -----

def main(args):

    # tbprofiler_results_location = 'tbprofiler_pakistan_results/'
    metadata_file = args.metadata_file
    id_key = args.id_key
    tbprofiler_results_location = args.tbp_results
    suffix = args.suffix
    outfile = args.outfile
    # db = args.db

    # -------------
    # READ IN DATA
    # -------------

    #  Read in metadata
    with open(metadata_file) as mf:
        metadata_reader = csv.DictReader(mf)
        meta_dict = {}
        for row in metadata_reader:
            # Make the id the key, but also recapitulate the id in the key-values by including everything
            meta_dict[row[id_key]] = row

    # json_files_list = ["{}{}".format(tbprofiler_results_location, i) for i in os.listdir(tbprofiler_results_location)]

    # mutations_dict = {}
    mutations_dict = defaultdict(dict)
    # for file in json_files_list:
    for samp in meta_dict:
        # Open the json file for the sample
        data = json.load(open("%s%s%s" % (tbprofiler_results_location, samp, suffix)))
        # samp = data["id"]
        inh_dst = meta_dict[samp]['isoniazid']
        lins = [lin['lin'] for lin in data['lineage']]
        lin = lins[len(lins) - 1] # I hate Python!
        for var in data["dr_variants"] + data["other_variants"]:
            # if "drugs" in var:
            #     drugs = [drug['drug'] for drug in var['drugs']]
            # else:
            #     drugs = "unknown"
            if var["gene"] == "ahpC" and var['type'] != 'synonymous' and "drugs" not in var and var["freq"] >= 0.7:

                mutations_dict[samp] = {'wgs_id': samp,
                'drtype':data["drtype"],
                'lineage':lin,
                'gene':var["gene"],
                'change':var["change"],
                'freq':var["freq"], 
                # 'drugs': drugs,   
                'inh_dst': inh_dst}

    fieldnames = tuple(next(iter(mutations_dict.values())).keys())

    with open(outfile, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        # Loop over the dictionaries, appending each dictionary as a row in the file
        for id in mutations_dict:
            writer.writerow(mutations_dict[id])

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata-file', default = '', type = str, help = 'metadata file name with column of sample ids - only processes samples in this file')
parser.add_argument('--id-key', default = '', type = str, help = 'column name in metadata file with sample ids')
# parser.add_argument('--db',default="tbdb",type=str,help='prefix to bed file for locus-DR associations')
parser.add_argument('--tbp-results', default="results/",type=str,help='tbprofiler results directory (json files)')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.add_argument('--outfile',default="ahpc_variants.txt",type=str,help='name of output file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)


