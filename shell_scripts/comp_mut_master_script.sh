#! /bin/bash

# comp_mut master script

set -e
set -u
set -o pipefail

# ----------------
# Setup - GLOBAL
# ----------------

cd ~/comp_mut


# ----------
# Variables
# ----------

# Directories

# Existing directories
tbp_results_dir=/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/
main_metadata_dir=../metadata/
local_metadata_dir=metadata/
local_json_dir=json/


# Files

# Existing files
main_tb_metadata_file=${main_metadata_dir}tb_data_18_02_2021.csv

# Created files
# head_metadata_file=${local_metadata_dir}head_metadata.csv # testing
test_sample_list=test_sample_list
temp_json_files_list=files.tmp



# Pull the samples from the metadata and only copy across json files for which there is metadata

# Get the sample IDs and put in a file
cat ${main_tb_metadata_file} | csvtk cut -f "wgs_id" | tail -n +2 | head -50 > ${test_sample_list}
# List the files in the TBprofiler results dir (with absolute path - find command) and grep the ones from the metadata
find ${tbp_results_dir} -type f | grep -f ${test_sample_list} > ${temp_json_files_list}
# Copy these files across 
cat ${temp_json_files_list} | parallel --bar -j1 cp {} ${local_json_dir}




# -----------------------------------------------------------
# Run tbprofiler_template.py to find novel mutations in katG
# -----------------------------------------------------------


# python python_scripts/tbprofiler_template_testing.py --sample_list_file test_sample_list --dir ../pakistan/tbprofiler_pakistan_results/json

python python_scripts/tbprofiler_template.py --dir ${tbp_results_dir}







