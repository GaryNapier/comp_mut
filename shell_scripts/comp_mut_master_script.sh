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

# Directories
# ------------

# Existing directories
tbp_results_dir=/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/
main_metadata_dir=../metadata/
local_metadata_dir=metadata/
database_dir=../pipeline/db/
results_dir=results/

# Files
# ------

# Existing files
main_tb_metadata_file=${main_metadata_dir}tb_data_18_02_2021.csv
known_comp_mut_file=${database_dir}compensatory_mutations.csv 

# Created files
head_metadata_file=${main_metadata_dir}head_metadata.csv # testing
sample_list=test_sample_list.txt
temp_json_files_list=files.tmp
# find_comp_mutations.py files
novel_comp_mut_data_file=${results_dir}novel_ahpc_mutations.txt
# filter_comp_mut.R files
model_results_file=${results_dir}ahpc_model_results.csv

# Variables 
drug_of_interest=isoniazid
id_key=wgs_id


# ------------------------
# Find all comp mutations
# ------------------------

if [ ! -f ${novel_comp_mut_data_file} ]; then
echo " --- PULLING COMPENSATORY MUTATIONS DATA FOR ${drug_of_interest} RUNNING python_scripts/find_novel_comp_mutations.py --- "
python python_scripts/find_novel_comp_mutations.py \
--metadata-file ${main_tb_metadata_file} \
--known-comp-mut-file ${known_comp_mut_file}
--drug-of-interest ${drug_of_interest} \
--id-key ${wgs_id} \
--tbp-results ${tbp_results_dir} \
--outfile ${novel_comp_mut_data_file}
fi

# ------------------------------------------------------------------------------------------------------------------------------
# Run filter_comp_mut.R - run GLM models to see if DST is significantly predicted against presence of each mutation and lineage
# ------------------------------------------------------------------------------------------------------------------------------

if [ ! -f ${model_results_file} ]; then
echo " --- FILTERING POTENTIAL NEW COMPENSATORY MUTATIONS - RUNNING r_scripts/filter_comp_mut.R --- "
Rscript r_scripts/filter_comp_mut.R \
--metadata_file ${main_tb_metadata_file} \
--drug_of_interest ${drug_of_interest}\
--novel_comp_mut_data_file ${novel_comp_mut_data_file} \
--outfile ${model_results_file} 
fi

# -----------------------------------------------------------

# -----------------------------------------------------------


