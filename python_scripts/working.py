import json
from collections import defaultdict, Counter
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
import tbprofiler
from csv import DictReader

# def main(args):

# mutations_file = args.mutations_file
# metadata_file = args.metadata_file
# id_key = args.id_key

mutations_file = "metadata/novel_ahpc_mutations.txt"
metadata_file = "../metadata/tb_data_18_02_2021.csv"
id_key = "wgs_id"
mutaions_key = "wgs_id" # column name of mutation names e.g. "c.-101A>G"

# -------------
# READ IN DATA
# -------------

# Read in mutations data
mutations_dict = {}
with open(mutations_file, 'r') as mf:
    mutations_reader = DictReader(mf, delimiter='\t')
    for row in mutations_reader:
        # Make the id the key, but also recapitulate the id in the key-values by including everything
        mutations_dict[row[mutaions_key]] = row

# Read in metadata
meta_dict = {}
with open(metadata_file) as mf:
    metadata_reader = csv.DictReader(mf)
    for row in metadata_reader:
        # Make the id the key, but also recapitulate the id in the key-values by including everything
        meta_dict[row[id_key]] = row


# ------------------
# Process mutations
# ------------------

# Get set of unique mutations
mutations_set = []
for key in mutations_dict:
    mutations_set.append(mutations_dict[key]['change'])
mutations_set = set(mutations_set)



parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata-file', default = '', type = str, help = 'metadata file name with column of sample ids - only processes samples in this file')
parser.add_argument('--mutations-file', default = '', type = str, help = 'txt file of novel mutations found with find_ahpc_mutations.py; columns: wgs_id, drtype, lineage, gene, change, freq, inh_dst)
parser.add_argument('--id-key', default = '', type = str, help = 'column name in metadata file with sample ids')
parser.add_argument('--outfile',default="ahpc_mutations_stats.txt",type=str,help='name of output file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)








# mutation 1
#         Has mutation?    INH_DST
# samp 1  0                   0
# samp 2  1                   1
# samp 3  0                   0


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# def get_gene2drugs(bed_file):
#     lt2drugs = {}
#     for l in open(bed_file):
#         row = l.strip().split()
#         lt2drugs[row[4]] = row[5].split(",")
#     return lt2drugs

# def get_conf_dict(library_prefix):
#     files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
#     conf = {}
#     for key in files:
#         sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
#         conf[key] = pp.filecheck(library_prefix+files[key])
#     return conf

# db = 'tbdb'

# conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % db)

# # Get a dictionary mapping the locus_tags to drugs: {"Rv1484": ["isoniazid","ethionamide"], ... etc. }
# # gene2drugs = tbprofiler.get_gene2drugs(conf["bed"])
# gene2drugs = get_gene2drugs(conf["bed"])

# metadata_file = "../metadata/tb_data_18_02_2021.csv"

# project_samples = ['SRR1723521']
# project_dir = "/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/"
# suffix = ".results.json"

# project_variants = defaultdict(list)
# # for s in tqdm(project_samples):
# s = project_samples[0]
# data = json.load(open(pp.filecheck("%s/%s%s" % (project_dir,s,suffix))))
# for var in data["dr_variants"] + data["other_variants"]:
#     if var["type"]=="synonymous": continue
#     if "Ser315Thr" in var["change"]:
#         print(var["change"])
#      # {('katG', 'p.Ser315Thr', 'isoniazid'): ['SRR1723521'],
#             #  ('gyrA', 'p.Glu21Gln', 'NA'): ['SRR1723521'],
#             #  ('gyrA', 'p.Ser95Thr', 'NA'): ['SRR1723521'],
#             #  ('gyrA', 'p.Gly668Asp', 'NA'): ['SRR1723521'],
#             #  ('mshA', 'p.Ala187Val', 'NA'): ['SRR1723521'],
#             #  ('rpsL', 'c.-165T>C', 'NA'): ['SRR1723521'],
#             #  ('katG', 'p.Arg463Leu', 'NA'): ['SRR1723521'],
#             #  ('ahpC', 'c.-101A>G', 'NA'): ['SRR1723521'],
#             #  ('thyA', 'p.Pro3Ala', 'NA'): ['SRR1723521'],
#             #  ('ald', 'c.-32T>C', 'NA'): ['SRR1723521'],
#             #  ('gid', 'p.Glu92Asp', 'NA'): ['SRR1723521']})
#     if "drugs" in var:
#         for d in var["drugs"]:
#             project_variants[(var["gene"],var["change"],d["drug"])].append(s)
#     else:
#         project_variants[(var["gene"],var["change"],"NA")].append(s)

# global_dir = project_dir
# global_samples = ['SRR1723521']
# global_variants = defaultdict(list)
# global_drtypes = defaultdict(list)
# global_variant_lineages = defaultdict(list)
# # for s in tqdm(global_samples):
# s = global_samples[0]
# data = json.load(open(pp.filecheck("%s/%s%s" % (global_dir,s,suffix))))
# for var in data["dr_variants"] + data["other_variants"]:
#     # defaultdict(list,
#     #         {('katG', 'p.Ser315Thr'): ['SRR1723521'],
#     #          ('gyrA', 'p.Glu21Gln'): ['SRR1723521'],
#     #          ('gyrA', 'p.Ser95Thr'): ['SRR1723521'],
#     #          ('gyrA', 'p.Gly668Asp'): ['SRR1723521'],
#     #          ('fgd1', 'c.960T>C'): ['SRR1723521'],
#     #          ('mshA', 'p.Ala187Val'): ['SRR1723521'],
#     #          ('rpoB', 'c.3225T>C'): ['SRR1723521'],
#     #          ('rpsL', 'c.-165T>C'): ['SRR1723521'],
#     #          ('rpsA', 'c.636A>C'): ['SRR1723521'],
#     #          ('tlyA', 'c.33A>G'): ['SRR1723521'],
#     #          ('katG', 'p.Arg463Leu'): ['SRR1723521'],
#     #          ('ahpC', 'c.-101A>G'): ['SRR1723521'],
#     #          ('thyA', 'c.666G>C'): ['SRR1723521'],
#     #          ('thyA', 'p.Pro3Ala'): ['SRR1723521'],
#     #          ('ald', 'c.-32T>C'): ['SRR1723521'],
#     #          ('embC', 'c.2781C>T'): ['SRR1723521'],
#     #          ('embA', 'c.228C>T'): ['SRR1723521'],
#     #          ('gid', 'c.615T>C'): ['SRR1723521'],
#     #          ('gid', 'p.Glu92Asp'): ['SRR1723521']})
#     global_variants[(var["gene"],var["change"])].append(s)
#     # defaultdict(list,
#     #         {('katG', 'p.Ser315Thr'): ['lineage2.2.1'],
#     #          ('gyrA', 'p.Glu21Gln'): ['lineage2.2.1'],
#     #          ('gyrA', 'p.Ser95Thr'): ['lineage2.2.1'],
#     #          ('gyrA', 'p.Gly668Asp'): ['lineage2.2.1'],
#     #          ('fgd1', 'c.960T>C'): ['lineage2.2.1'],
#     #          ('mshA', 'p.Ala187Val'): ['lineage2.2.1'],
#     #          ('rpoB', 'c.3225T>C'): ['lineage2.2.1'],
#     #          ('rpsL', 'c.-165T>C'): ['lineage2.2.1'],
#     #          ('rpsA', 'c.636A>C'): ['lineage2.2.1'],
#     #          ('tlyA', 'c.33A>G'): ['lineage2.2.1'],
#     #          ('katG', 'p.Arg463Leu'): ['lineage2.2.1'],
#     #          ('ahpC', 'c.-101A>G'): ['lineage2.2.1'],
#     #          ('thyA', 'c.666G>C'): ['lineage2.2.1'],
#     #          ('thyA', 'p.Pro3Ala'): ['lineage2.2.1'],
#     #          ('ald', 'c.-32T>C'): ['lineage2.2.1'],
#     #          ('embC', 'c.2781C>T'): ['lineage2.2.1'],
#     #          ('embA', 'c.228C>T'): ['lineage2.2.1'],
#     #          ('gid', 'c.615T>C'): ['lineage2.2.1'],
#     #          ('gid', 'p.Glu92Asp'): ['lineage2.2.1']})
#     global_variant_lineages[(var["gene"],var["change"])].append(data["sublin"])
# # defaultdict(list, {'Pre-MDR': ['SRR1723521']})
# global_drtypes[data["drtype"]].append(s)

# meta = {}
# for row in csv.DictReader(open(metadata_file)):
#     meta[row["wgs_id"]] = row

# # ---


# # def calculate_stats(gene,change,drug,dr_pos,project_freq):


# # for drtype in ["Sensitive","Pre-MDR","MDR","Pre-XDR","XDR","Other"]:
# #     samps = global_drtypes[drtype]
# #     res[drtype+"_freq"] = sum([1 for s in samps if s in g_samples])/len(samps)
# # return res


# results = []
# # for gene,change,drug in tqdm(project_variants):

# # gene = 'katG'
# # change = 'p.Ser315Thr'
# # drug = 'isoniazid'

# gene = 'gid'
# change = 'p.Glu92Asp' 
# drug = 'NA'


# project_freq = len(project_variants[(gene,change,drug)])

# if gene not in gene2drugs:
#     results.append({"gene":gene, 
#         "change":change, 
#         "drug":"NA", 
#         "dr_pos":"No",
#         "project_freq":project_freq,
#         "global_freq":"NA",
#         "drug_resistant_%":"NA"})


# # if drug=="NA": 
# possible_drugs = gene2drugs[gene]
# # for d in possible_drugs:

# # ---------------------------
# # calculate_stats() function
# # ---------------------------
# # def   calculate_stats(gene,change,drug,dr_pos,project_freq):
# # var = calculate_stats(gene,change,d   ,"No"  ,project_freq)
# d = possible_drugs[0]
# dr_pos = "No"

# g_samples = set(global_variants[(gene,change)])
# global_dsts = [int(meta[s][d]) for s in g_samples if meta[s][d]!="NA"] if d in list(list(meta.values())[0]) else []
# res = {
#     "gene":gene,
#     "change":change,
#     "drug":drug,
#     "dr_pos":dr_pos,
#     "project_freq":project_freq,
#     "global_N":len(g_samples),
#     "drug_resistant_%":"%.2f" % (sum(global_dsts)/len(global_dsts)) if len(global_dsts)!=0 else "NA"
# }



# for drtype in ["Sensitive","Pre-MDR","MDR","Pre-XDR","XDR","Other"]:
#     samps = global_drtypes[drtype]
#     res[drtype+"_freq"] = sum([1 for s in samps if s in g_samples])/len(samps) if len(samps) > 0 else 0

#     var["lineages"] = resolve_lineages(Counter(global_variant_lineages[(gene,change)]))
#     var["num_lineages"] = len(var["lineages"])
#     results.append(var)


# # else:
# #     var = calculate_stats(gene,change,drug,"Yes",project_freq)
# #     var["lineages"] = resolve_lineages(Counter(global_variant_lineages[(gene,change)]))
# #     var["num_lineages"] = len(var["lineages"])
# #     results.append(var)



