#! /bin/bash

# comp_mut master script

set -e
set -u
set -o pipefail

# ------
# Setup 
# ------

cd ~/comp_mut

# ----------
# Variables
# ----------

drug_of_interest=isoniazid
id_key=wgs_id

# ------------
# Directories
# ------------

# Existing directories
tbp_results_dir=/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/
main_metadata_dir=../metadata/
local_metadata_dir=metadata/
database_dir=../pipeline/db/
tbdb_dir=../tbdb/
results_dir=results/

# ------
# Files
# ------

# Existing files
main_tb_metadata_file=${main_metadata_dir}tb_data_18_02_2021.csv
known_comp_mut_file=${database_dir}compensatory_mutations.csv 
tbdb_file=${tbdb_dir}tbdb.csv
dr_types_file=${database_dir}dr_types.json
vars_exclude_file=${local_metadata_dir}var_exclude_katg_comp_mut.csv

# Created files
head_metadata_file=${main_metadata_dir}head_metadata.csv # testing
sample_list=test_sample_list.txt
temp_json_files_list=files.tmp
# find_comp_mutations.py files
novel_comp_mut_data_file=${results_dir}${drug_of_interest}_novel_comp_mut_data.txt
# filter_comp_mut.R files
novel_comp_mut_model_results_file=${results_dir}${drug_of_interest}_novel_comp_mut_model_results.csv
# comp_mut2res_mut.py files
potential_res_mutations_outfile=${results_dir}${drug_of_interest}_potential_res_mutations.csv


# ---------------------------------------------------
# Find all novel comp mutations for drug of interest
# ---------------------------------------------------

if [ ! -f ${novel_comp_mut_data_file} ]; then
echo " --- PULLING COMPENSATORY MUTATIONS DATA FOR ${drug_of_interest} RUNNING python_scripts/find_novel_comp_mutations.py --- "
python python_scripts/find_novel_comp_mutations.py \
--metadata-file ${main_tb_metadata_file} \
--known-comp-mut-file ${known_comp_mut_file} \
--drug-of-interest ${drug_of_interest} \
--id-key ${id_key} \
--tbp-results ${tbp_results_dir} \
--outfile ${novel_comp_mut_data_file}
fi

# ------------------------------------------------------------------------------------------------------------------------------------
# Run filter_novel_comp_mut.R - run GLM models to see if DST is significantly predicted against presence of each mutation and lineage
# ------------------------------------------------------------------------------------------------------------------------------------

if [ ! -f ${novel_comp_mut_model_results_file} ]; then
echo " --- FILTERING POTENTIAL NEW COMPENSATORY MUTATIONS - RUNNING r_scripts/filter_novel_comp_mut.R --- "
Rscript r_scripts/filter_novel_comp_mut.R \
--metadata_file ${main_tb_metadata_file} \
--novel_comp_mut_data_file ${novel_comp_mut_data_file} \
--drug_of_interest ${drug_of_interest} \
--outfile ${novel_comp_mut_model_results_file} 
fi

# -----------------------------------------------------------
# Filter potential novel compensatory mutations with tbprofiler critera (n lineages, proportion of samples drug resistant, proportion of samples DST resistant)
# 
# -----------------------------------------------------------

python python/comp_mut2res_mut.py \

--potential-comp-mut-file ${novel_comp_mut_model_results_file} \
--metadata-file ${main_tb_metadata_file} \
--tbdb-file ${tbdb_file} \
--drtypes-file ${dr_types_file} \ 
--known-comp-mut-file ${known_comp_mut_file} \
--tbprofiler-results-dir ${tbp_results_dir} \
--vars-exclude-file ${vars_exclude_file} \
--potential-res-mut-outfile ${potential_res_mutations_outfile} \
--drug-of-interest ${drug_of_interest}

# --comp-mut-genes 

