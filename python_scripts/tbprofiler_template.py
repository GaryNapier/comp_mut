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
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf

def main(args):
    # Get a dictionary with the database file: {"ref": "/path/to/fasta" ... etc. }
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)

    # Get a dictionary mapping the locus_tags to drugs: {"Rv1484": ["isoniazid","ethionamide"], ... etc. }
    locus_tag2drugs = tbprofiler.get_lt2drugs(conf["bed"])

    # If a list of samples is supplied through the args object, store it in a list else get the list from looking in the results direcotry
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    # Loop through the sample result files
    genes_to_ananalyse = set(["ahpC","katG"])
    sample_data = defaultdict(dict)
    variant2samples = defaultdict(set)
    potential_other_mutations  = []
    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        sample_data[s]["drtype"] = data["drtype"]
        dr_variants = []
        other_variants = []
        for var in data["dr_variants"] + data["other_variants"]:
            v = (var["gene"],var["change"])
            if var["gene"] in genes_to_ananalyse:
                variant2samples[v].add(s)
                if "drugs" in var:
                    dr_variants.append(v)
                else:
                    other_variants.append(v)
        
        dr_variant_genes = [v[0] for v in dr_variants]
        
        if "ahpC" in dr_variant_genes and "katG" not in dr_variant_genes:
            print(s)
            for var in other_variants:
                potential_other_mutations.append(var)

    
    potential_other_mutations = Counter(potential_other_mutations)
    for gene,change in potential_other_mutations:
        samples_with_mut = variant2samples[(gene,change)]
        drtypes = [sample_data[s]["drtype"] for s in samples_with_mut]
        print(gene,change,potential_other_mutations[(gene,change)],Counter(drtypes))



# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
