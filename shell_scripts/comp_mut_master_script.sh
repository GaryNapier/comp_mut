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
all_ahpc_mutations_file=${local_metadata_dir}all_ahpc_mutations.txt

# ------------------------
# Find all ahpC mutations
# ------------------------

echo "RUNNING python_scripts/find_ahpc_mutations.py"
python python_scripts/find_ahpc_mutations.py --metadata-file ${main_tb_metadata_file} --id-key wgs_id --tbp-results ${tbp_results_dir} --outfile ${all_ahpc_mutations_file}

# -----------------------------------------------------------
# Run tbprofiler_template.py to find novel mutations in katG
# -----------------------------------------------------------


# python python_scripts/tbprofiler_template_testing.py --sample_list_file sample_list --dir ../pakistan/tbprofiler_pakistan_results/json

# python python_scripts/tbprofiler_template.py --dir ${tbp_results_dir}

