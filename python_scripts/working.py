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
