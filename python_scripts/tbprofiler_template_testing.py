#! /usr/bin/env python

# Load useful libraries
import json
from collections import defaultdict, Counter
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
import tbprofiler


def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        # sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        print("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf

db = "tbdb"
conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % db)

locus_tag2drugs = tbprofiler.get_lt2drugs(conf["bed"])

# If a list of samples is supplied through the args object, store it in a list else get the list from looking in the results direcotry

# samples = "test_sample_list"
# args_dir = "../pakistan/tbprofiler_pakistan_results/json"
args_dir = "~jody/tbprofiler_tests/results/"
suffix = ".results.json"

# if args.samples:
if samples:
    samples = [x.rstrip() for x in open(samples).readlines()]
else:
    samples = [x.replace(suffix,"") for x in os.listdir(args_dir) if x[-len(suffix):]==suffix]


# Loop through the sample result files
genes_to_analyse = set(["ahpC","katG"])
sample_data = defaultdict(dict)
variant2samples = defaultdict(set)
potential_other_mutations  = []

# for s in samples:
s = samples[0]
# Data has the same structure as the .result.json files
# Open the json file
data = json.load(open("%s/%s%s" % (args_dir,s,suffix)))
# Get the DR type ('MDR', 'XDR', etc) and store as 'drtype':<type> dict as the value of the samples name key:
# {'ERR688030': {'drtype': 'XDR'},
#  'ERR688043': {'drtype': 'XDR'},
#  'ERR1213917': {'drtype': 'Sensitive'}}
sample_data[s]["drtype"] = data["drtype"]

dr_variants = []
other_variants = []
# Loop through the concatenated list of DR variants and other variants, pull out the gene and change, 
# store as tuple
# for var in data["dr_variants"] + data["other_variants"]:

dr_other = data["dr_variants"] + data["other_variants"]
var = [x for x in dr_other if x["gene"] in genes_to_analyse]
var = var[0]

v = (var["gene"], var["change"])
# Add to dict if one of the genes in the list above

# if var["gene"] in genes_to_analyse:

# Store as dict - <gene/change tuple> : <sample[type = set]>
# e.g. 
# {('katG', 'p.Ser315Thr'): {'ERR2510563'},
# ('katG', 'p.Arg463Leu'): {'ERR2510563'},
# ('ahpC', 'c.-88G>A'): {'ERR2510563'}}
variant2samples[v].add(s)
# Separate the known DR variants and 'other' and save the gene/change tuple
if "drugs" in var:
    # e.g. : [('katG', 'p.Ser315Thr')]
    dr_variants.append(v)
else:
    # e.g. : [('katG', 'p.Arg463Leu'), ('ahpC', 'c.-88G>A')]
    other_variants.append(v)

# ['katG']
dr_variant_genes = [gene[0] for gene in dr_variants]

if "ahpC" in dr_variant_genes and "katG" not in dr_variant_genes:
    print(s)
    for variant in other_variants:
        print(variant)
        # [('katG', 'p.Arg463Leu'), ('ahpC', 'c.-88G>A')]
        potential_other_mutations.append(varariant)

potential_other_mutations = Counter(potential_other_mutations)
for gene,change in potential_other_mutations:
    samples_with_mut = variant2samples[(gene,change)]
    drtypes = [sample_data[s]["drtype"] for s in samples_with_mut]
    print(gene,change,potential_other_mutations[(gene,change)],Counter(drtypes))





