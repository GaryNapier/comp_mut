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

mf.close()

meta_dict

# SAMEA2534433.results.json
# ERR1465765.results.json
# SAMEA3715554.results.json
# ERR1465906.results.json
# ERR1193661.results.json
# ERR1193674.results.json
# ERR1193709.results.json
# ERR1465895.results.json
# ERR757145.results.json

files_dir = "json/"

json_files_list = ["{}{}".format(files_dir, i) for i in os.listdir(files_dir)]

# mutations_dict = {}
mutations_dict = defaultdict(dict)
for file in json_files_list:
    # Open the json file for the sample
    data = json.load(open(file))
    for var in data["dr_variants"] + data["other_variants"]:
        if var["gene"] == "ahpC":
            samp = data["id"]
            mutations_dict[samp]["drtype"] = data["drtype"]
            mutations_dict[samp]["lineage"] = data["lineage"]
            mutations_dict[samp]["gene"] = var["gene"]
            mutations_dict[samp]["change"] = var["change"]



