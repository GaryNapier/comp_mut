import json
from collections import defaultdict, Counter
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
import tbprofiler

# metadata = "../metadata/tb_data_18_02_2021.csv"
metadata_file = "metadata/head_metadata.csv"

# Variables
id_key = "wgs_id"

#  Read in metadata
with open(metadata_file) as mf:
    metadata_reader = csv.DictReader(mf)
    meta_dict = {}
    for row in metadata_reader:
        # Make the id the key, but also recapitulate the id in the key-values by including everything
        meta_dict[row[id_key]] = row


# files_dir = "json/"

# json_files_list = ["{}{}".format(files_dir, i) for i in os.listdir(files_dir)]

# json_files_list = ["json/SRR11662391.results.json"]

tbprofiler_results_location = "/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/"
suffix = ".results.json"


# mutations_dict = {}
mutations_dict = defaultdict(dict)
# for file in json_files_list:
for samp in meta_dict:
    # Open the json file for the sample
    # data = json.load(open(file))
    data = json.load(open("%s%s%s" % (tbprofiler_results_location, samp, suffix)))
    # samp = data["id"]
    inh_dst = meta_dict[samp]['isoniazid']
    lins = [lin['lin'] for lin in data['lineage']]
    lin = lins[len(lins) - 1] # I hate Python!
    for var in data["dr_variants"] + data["other_variants"]:
        if "drugs" in var:
            drugs = [drug['drug'] for drug in var['drugs']]
        else:
            drugs = "unknown"
        if var["gene"] == "ahpC" and var['type'] != 'synonymous':
            mutations_dict[samp] = {'wgs_id': samp,
            'drtype':data["drtype"],
            'lineage':lin,
            'gene':var["gene"],
            'change':var["change"],
            'freq':var["freq"], 
            'drugs': drugs, 
            'inh_dst': inh_dst}


fieldnames = tuple(next(iter(mutations_dict.values())).keys())

outfile = 'mut_test.txt'
with open(outfile, 'w') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    # Loop over the dictionaries, appending each dictionary as a row in the file
    for id in mutations_dict:
        writer.writerow(mutations_dict[id])