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

# Existing directories
tbp_results_dir=/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/
main_metadata_dir=../metadata/
local_metadata_dir=metadata/

# Files

# Existing files
main_tb_metadata_file=${main_metadata_dir}tb_data_18_02_2021.csv

# Created files

head_metadata_file=${main_metadata_dir}head_metadata.csv # testing
sample_list=test_sample_list.txt
temp_json_files_list=files.tmp
# find_ahpc_mutations.py
novel_ahpc_mutations_file=${local_metadata_dir}novel_ahpc_mutations.txt
# ahpc_mut_glm.R
ahpc_model_results_file=${local_metadata_dir}ahpc_model_results.csv


# ------------------------
# Find all ahpC mutations
# ------------------------

if [ ! -f ${novel_ahpc_mutations_file} ]; then
echo " --- RUNNING python_scripts/find_ahpc_mutations.py --- "
python python_scripts/find_ahpc_mutations.py --metadata-file ${main_tb_metadata_file} --id-key wgs_id --tbp-results ${tbp_results_dir} --outfile ${novel_ahpc_mutations_file}
fi

# ------------------------------------------------------------------------------------------------------------
# Run ahpc_mut_glm.R to run models to see if DST is significantly predicted against each mutation and lineage
# ------------------------------------------------------------------------------------------------------------

if [ ! -f ${ahpc_model_results_file} ]; then
echo " --- RUNNING r_scripts/ahpc_mut_glm.R --- "
Rscript r_scripts/ahpc_mut_glm.R --metadata_file ${main_tb_metadata_file} --mutations_data_file ${novel_ahpc_mutations_file} --outfile ${ahpc_model_results_file}
fi

# -----------------------------------------------------------
# Run tbprofiler_template.py to find novel mutations in katG
# -----------------------------------------------------------


# python python_scripts/tbprofiler_template_testing.py --sample_list_file sample_list --dir ../pakistan/tbprofiler_pakistan_results/json

# python python_scripts/tbprofiler_template.py --dir ${tbp_results_dir}

