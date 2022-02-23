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
from collections import Counter
import requests
from contextlib import closing
import re
from python_scripts.utils import *

# 'SAMEA2534433'

ahpc_glm_results_file = "metadata/ahpc_model_results.csv"
metadata_file = "../metadata/tb_data_18_02_2021.csv"
tbdb_file = "../tbdb/tbdb.csv"
drtypes_file = "../pipeline/db/dr_types.json"
# tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/'
# tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/gatk/results'
tbprofiler_results_dir = '/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/'
# metadata_id_key = "wgs_id"
suffix = ".results.json"
genes = ('ahpC', 'katG', 'fabG1')
vars_exclude_file = 'metadata/var_exclude_katg_comp_mut.csv'
drug_of_interest = 'isoniazid'


# TESTING
# file = "%s/%s%s" % (tbprofiler_results_dir, 'SAMEA2534433', suffix)
# data = json.load(open(file))

# ------------------------------------------------------------
# KATG - Co-occurrence with fabG1 and comparison to Ser315Thr
# ------------------------------------------------------------

# potential_res_mut_filtered
# mutation2sample
# sample2mutation
# meta_dict
# drug_of_interest

co_gene = 'fabG1'

co_gene_variants = set()
for var in potential_res_mut_filtered:
    samps = mutation2sample[var]
    
    for samp in samps:
        variants = sample2mutation[samp]
        co_gene_variants = [v for v in variants if v[0]==co_gene]
        print(co_gene_variants)


# [var for var in mutation2sample if var[0] == co_gene]

# co_gene_vars = [var for var in mutation2sample if var[0] == co_gene]

# co_gene_samps = [mutation2sample[var] for var in co_gene_vars]

# Take all_katg - aggregate samples - how many have fabG1?


# Compare to p.Ser315Thr - get proportion of samples with p.Ser315Thr that don't have a fabG1


Ser315Thr_fabg_samps = []
Ser315Thr_no_fabg_samps = []
for samp in all_data:
    mutations = all_data[samp]['mutations']
    if any((x['gene'] == 'katG' and x['change'] == 'p.Ser315Thr') and any(x['gene'] == 'fabG1' for x in mutations) for x in mutations):
        Ser315Thr_fabg_samps.append(samp)
    if any((x['gene'] == 'katG' and x['change'] == 'p.Ser315Thr') and not any(x['gene'] == 'fabG1' for x in mutations) for x in mutations):
        Ser315Thr_no_fabg_samps.append(samp)

Ser315Thr_fabg_prop = round(len(Ser315Thr_fabg_samps) / (len(Ser315Thr_fabg_samps) + len(Ser315Thr_no_fabg_samps)), 3)
# 0.158


# ================================================================================


# find_novel_comp_mutations.py edits

# Need to pull comp mutations data from jsons using the drug of interest as input
# Take the comp mutations file, get the genes by filtering the drug of interest col 



# ================================================================================


# Need to get all the genes associated with isoniazid

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
from collections import Counter
import requests
from contextlib import closing
import re
from python_scripts.utils import *

drug_of_interest = 'isoniazid'


tbdb_file = "../tbdb/tbdb.csv"
# Read in tbdb file
with open(tbdb_file, 'r') as f:
    tbdb_dict = csv_to_dict_multi(f, 'Drug')


for drug in tbdb_dict:
    print(tbdb_dict[drug])


set([l.strip().split()[0] for l in open(genes_file)])







fst_results_url = 'https://raw.githubusercontent.com/GaryNapier/tb-lineages/main/fst_results_clean_fst_1_for_paper.csv'
# See https://www.codegrepper.com/code-examples/python/how+to+read+a+csv+file+from+a+url+with+python for pulling data from url
with closing(requests.get(fst_results_url, stream=True)) as r:
    f = (line.decode('utf-8') for line in r.iter_lines())
    fst_dict = csv_to_dict_multi(f)

lin_specific_variants = []
for gene in fst_dict:
    if gene in genes:
        for var in fst_dict[gene]:
            lin_specific_variants.append( tuple( [ gene, reformat_mutations(var['aa_pos']) ] ) )

# Read in variants to be excluded
vars_exclude = []
for l in open(vars_exclude_file):
    vars_exclude.append(tuple(l.strip().split(',')))

# Concat
vars_exclude = vars_exclude + lin_specific_variants 