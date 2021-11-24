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

    # Get the gene, change, drug combinations for each sample.
    # Example:
    # {('katG', 'p.Ser315Thr', 'isoniazid'): ['SRR1723521'],
    #  ('gyrA', 'p.Glu21Gln', 'NA'): ['SRR1723521'],
    #  ('gyrA', 'p.Ser95Thr', 'NA'): ['SRR1723521'],
    #  ('gyrA', 'p.Gly668Asp', 'NA'): ['SRR1723521'],
    #  ('mshA', 'p.Ala187Val', 'NA'): ['SRR1723521'],
    #  ('rpsL', 'c.-165T>C', 'NA'): ['SRR1723521'],
    #  ('katG', 'p.Arg463Leu', 'NA'): ['SRR1723521'],
    #  ('ahpC', 'c.-101A>G', 'NA'): ['SRR1723521'],
    #  ('thyA', 'p.Pro3Ala', 'NA'): ['SRR1723521'],
    #  ('ald', 'c.-32T>C', 'NA'): ['SRR1723521'],
    #  ('gid', 'p.Glu92Asp', 'NA'): ['SRR1723521']})
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

    # Same for global samples, except no filter for synonymous, no 'drug' and also get the gene-change combination for lineage:
    # defaultdict(list,
    #         {('katG', 'p.Ser315Thr'): ['SRR1723521'],
    #          ('gyrA', 'p.Glu21Gln'): ['SRR1723521'],
    #          ('gyrA', 'p.Ser95Thr'): ['SRR1723521'],
    #          ('gyrA', 'p.Gly668Asp'): ['SRR1723521'],
    #          ('fgd1', 'c.960T>C'): ['SRR1723521'],
    #          ('mshA', 'p.Ala187Val'): ['SRR1723521'],
    #          ('rpoB', 'c.3225T>C'): ['SRR1723521'],
    #          ('rpsL', 'c.-165T>C'): ['SRR1723521'],
    #          ('rpsA', 'c.636A>C'): ['SRR1723521'],
    #          ('tlyA', 'c.33A>G'): ['SRR1723521'],
    #          ('katG', 'p.Arg463Leu'): ['SRR1723521'],
    #          ('ahpC', 'c.-101A>G'): ['SRR1723521'],
    #          ('thyA', 'c.666G>C'): ['SRR1723521'],
    #          ('thyA', 'p.Pro3Ala'): ['SRR1723521'],
    #          ('ald', 'c.-32T>C'): ['SRR1723521'],
    #          ('embC', 'c.2781C>T'): ['SRR1723521'],
    #          ('embA', 'c.228C>T'): ['SRR1723521'],
    #          ('gid', 'c.615T>C'): ['SRR1723521'],
    #          ('gid', 'p.Glu92Asp'): ['SRR1723521']})

    # and

        # defaultdict(list,
    #         {('katG', 'p.Ser315Thr'): ['lineage2.2.1'],
    #          ('gyrA', 'p.Glu21Gln'): ['lineage2.2.1'],
    #          ('gyrA', 'p.Ser95Thr'): ['lineage2.2.1'],
    #          ('gyrA', 'p.Gly668Asp'): ['lineage2.2.1'],
    #          ('fgd1', 'c.960T>C'): ['lineage2.2.1'],
    #          ('mshA', 'p.Ala187Val'): ['lineage2.2.1'],
    #          ('rpoB', 'c.3225T>C'): ['lineage2.2.1'],
    #          ('rpsL', 'c.-165T>C'): ['lineage2.2.1'],
    #          ('rpsA', 'c.636A>C'): ['lineage2.2.1'],
    #          ('tlyA', 'c.33A>G'): ['lineage2.2.1'],
    #          ('katG', 'p.Arg463Leu'): ['lineage2.2.1'],
    #          ('ahpC', 'c.-101A>G'): ['lineage2.2.1'],
    #          ('thyA', 'c.666G>C'): ['lineage2.2.1'],
    #          ('thyA', 'p.Pro3Ala'): ['lineage2.2.1'],
    #          ('ald', 'c.-32T>C'): ['lineage2.2.1'],
    #          ('embC', 'c.2781C>T'): ['lineage2.2.1'],
    #          ('embA', 'c.228C>T'): ['lineage2.2.1'],
    #          ('gid', 'c.615T>C'): ['lineage2.2.1'],
    #          ('gid', 'p.Glu92Asp'): ['lineage2.2.1']})

    global_variants = defaultdict(list)
    global_drtypes = defaultdict(list)
    global_variant_lineages = defaultdict(list)
    for s in tqdm(global_samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.global_dir,s,args.suffix))))
        for var in data["dr_variants"] + data["other_variants"]:
            global_variants[(var["gene"],var["change"])].append(s)
            global_variant_lineages[(var["gene"],var["change"])].append(data["sublin"])

        # Store the DR type for each samp:
        # Example:
        # defaultdict(list, {'Pre-MDR': ['SRR1723521']})
        global_drtypes[data["drtype"]].append(s)
    
    meta = {}
    for row in csv.DictReader(open(args.meta)):
        meta[row["wgs_id"]] = row

    def calculate_stats(gene,change,drug,dr_pos,project_freq):

        # Get the samples that have the combination of gene & change from the global variants list
        g_samples = set(global_variants[(gene,change)])
        # Extract the DST status (0/1) as numeric values, ignoring missing values
        # Final DST freq will be proportion of sensitive and susceptible for just those samples that have been tested. 
        global_dsts = [int(meta[s][drug]) for s in g_samples if meta[s][drug]!="NA"] if drug in list(list(meta.values())[0]) else []
        # Store in dict
        res = {
            "gene":gene,"change":change,"drug":drug,"dr_pos":dr_pos,
            "project_freq":project_freq,
            "global_N":len(g_samples),
            "drug_resistant_%":"%.2f" % (sum(global_dsts)/len(global_dsts)) if len(global_dsts)!=0 else "NA"
        }
        # Loop over DR types and get the number of samples in each category
        for drtype in ["Sensitive","Pre-MDR","MDR","Pre-XDR","XDR","Other"]:
            samps = global_drtypes[drtype]
            res[drtype+"_freq"] = sum([1 for s in samps if s in g_samples])/len(samps)
        return res

    results = []
    for gene,change,drug in tqdm(project_variants):
        # Get the freq for each gene,change,drug in the project (subset the dict on the gene,change,drug and get the number of samples (len))
        project_freq = len(project_variants[(gene,change,drug)])

        # If the gene is not a drug resistance-associated one, then make this the results dict entry:
        if gene not in gene2drugs:
            results.append({"gene":gene, 
            "change":change, 
            "drug":"NA",
            "dr_pos":"No",
            "project_freq":project_freq,
            "global_freq":"NA",
            "drug_resistant_%":"NA"})
        
        # If the drug is "NA" (the CHANGE is non-DR-associated), then... 
        elif drug=="NA":
            # look up the gene and get the drug(s)
            possible_drugs = gene2drugs[gene]
            # Loop over the drug(s) and 
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
# parser.add_argument('--project-samples',type=str,help='File with samples')
parser.add_argument('--global-samples',type=str,help='File with samples')
parser.add_argument('--out',type=str,help='directory for results',required=True)
parser.add_argument('--meta',type=str,help='metadata file',required=True)
# parser.add_argument('--project-dir',type=str,help='Directory containing results',required=True)
parser.add_argument('--global-dir',type=str,help='location of tbprofiler json files',required=True)
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
