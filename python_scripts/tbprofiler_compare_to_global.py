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

def resolve_lineages(data):
    data = dict(data)
    while True:
        improved = False
        for child in sorted(data,reverse= True,key=lambda x : len(x)):
            parent = ".".join(child.split(".")[:-1])
            if any([parent == x for x in data]):
                data[parent] = data[parent] + data[child] 
                del data[child]
                improved = True
                break
        if not improved:
            break
    return data

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
    gene2drugs = tbprofiler.get_gene2drugs(conf["bed"])

    # If a list of samples is supplied through the args object, store it in a list else get the list from looking in the results direcotry
    if args.project_samples:
        project_samples = [x.rstrip() for x in open(args.project_samples).readlines()]
    else:
        project_samples = [x.replace(args.suffix,"") for x in os.listdir(args.project_dir) if x[-len(args.suffix):]==args.suffix]

    if args.global_samples:
        global_samples = [x.rstrip() for x in open(args.global_samples).readlines()]
    else:
        global_samples = [x.replace(args.suffix,"") for x in os.listdir(args.global_dir) if x[-len(args.suffix):]==args.suffix]

    project_variants = defaultdict(list)
    for s in tqdm(project_samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.project_dir,s,args.suffix))))
        for var in data["dr_variants"] + data["other_variants"]:
            if var["type"]=="synonymous": continue
            if "Ser315Thr" in var["change"]:
                print(var["change"])
            if "drugs" in var:
                for d in var["drugs"]:
                    project_variants[(var["gene"],var["change"],d["drug"])].append(s)
            else:
                project_variants[(var["gene"],var["change"],"NA")].append(s)

    global_variants = defaultdict(list)
    global_drtypes = defaultdict(list)
    global_variant_lineages = defaultdict(list)
    for s in tqdm(global_samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.global_dir,s,args.suffix))))
        for var in data["dr_variants"] + data["other_variants"]:
            global_variants[(var["gene"],var["change"])].append(s)
            global_variant_lineages[(var["gene"],var["change"])].append(data["sublin"])
        global_drtypes[data["drtype"]].append(s)
    meta = {}
    for row in csv.DictReader(open(args.meta)):
        meta[row["wgs_id"]] = row

    def calculate_stats(gene,change,drug,dr_pos,project_freq):

        g_samples = set(global_variants[(gene,change)])
        global_dsts = [int(meta[s][drug]) for s in g_samples if meta[s][drug]!="NA"] if drug in list(list(meta.values())[0]) else []
        res = {
            "gene":gene,"change":change,"drug":drug,"dr_pos":dr_pos,
            "project_freq":project_freq,
            "global_freq":len(g_samples),
            "drug_resistant_%":"%.2f" % (sum(global_dsts)/len(global_dsts)) if len(global_dsts)!=0 else "NA"
        }
        for drtype in ["Sensitive","Pre-MDR","MDR","Pre-XDR","XDR","Other"]:
            samps = global_drtypes[drtype]
            res[drtype+"_freq"] = sum([1 for s in samps if s in g_samples])/len(samps)
        return res

    results = []
    for gene,change,drug in tqdm(project_variants):
        project_freq = len(project_variants[(gene,change,drug)])
        if gene not in gene2drugs:
            results.append({"gene":gene,"change":change,"drug":"NA","dr_pos":"No","project_freq":project_freq,"global_freq":"NA","drug_resistant_%":"NA"})
        elif drug=="NA":
            possible_drugs = gene2drugs[gene]
            for d in possible_drugs:
                var = calculate_stats(gene,change,d,"No",project_freq)
                var["lineages"] = resolve_lineages(Counter(global_variant_lineages[(gene,change)]))
                var["num_lineages"] = len(var["lineages"])
                results.append(var)
        else:
            var = calculate_stats(gene,change,drug,"Yes",project_freq)
            var["lineages"] = resolve_lineages(Counter(global_variant_lineages[(gene,change)]))
            var["num_lineages"] = len(var["lineages"])
            results.append(var)

    with open(args.out,"w") as O:
        writer = csv.DictWriter(O,fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--project-samples',type=str,help='File with samples')
parser.add_argument('--global-samples',type=str,help='File with samples')
parser.add_argument('--out',type=str,help='Directory containing results',required=True)
parser.add_argument('--meta',type=str,help='Directory containing results',required=True)
parser.add_argument('--project-dir',type=str,help='Directory containing results',required=True)
parser.add_argument('--global-dir',type=str,help='Directory containing results',required=True)
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
