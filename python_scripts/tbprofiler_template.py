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
    # Dictionary of file types and their suffixes
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
    # ahpC	isoniazid
    # katG	isoniazid
    genes_to_analyse = set(["ahpC","katG"])
    sample_data = defaultdict(dict)
    variant2samples = defaultdict(set)
    potential_other_mutations  = []
    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        # Open the json file for the sample
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        # Get the DR type ('MDR', 'XDR', etc) and store as 'drtype':<type> dict as the value of the samples name key:
        # {'ERR688030': {'drtype': 'XDR'},
        #  'ERR688043': {'drtype': 'XDR'},
        #  'ERR1213917': {'drtype': 'Sensitive'}}
        sample_data[s]["drtype"] = data["drtype"]
        
        dr_variants = []
        other_variants = []
        # Loop through the concatenated list of DR variants and other variants, pull out the gene and mutation change, 
        # store as tuple
        for var in data["dr_variants"] + data["other_variants"]:
            v = (var["gene"],var["change"])
            if var["gene"] in genes_to_analyse:
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
                # e.g. [('katG', 'p.Arg463Leu'), ('ahpC', 'c.-88G>A')]
                potential_other_mutations.append(variant)

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


# tbprofiler_template.py --samples <samples>