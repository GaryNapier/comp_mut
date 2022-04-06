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
vcf_remote=/mnt/storage7/jody/tb_ena/per_sample/
main_metadata_dir=../metadata/
local_metadata_dir=metadata/
database_dir=../pipeline/db/
tbdb_dir=../tbdb/
results_dir=results/
fasta_dir=fasta/
newick_dir=${results_dir}newick/
vcf_db=~/vcf/
vcf_dir=vcf/

# ------
# Files
# ------

# Existing files
main_tb_metadata_file=${main_metadata_dir}tb_data_18_02_2021.csv
KCM_file=${database_dir}compensatory_mutations.csv 
tbdb_file=${tbdb_dir}tbdb.csv
dr_types_file=${database_dir}dr_types.json
vars_exclude_file=${local_metadata_dir}var_exclude_comp_mut.csv
tc_file=${local_metadata_dir}${drug_of_interest}_tc.txt

# Created files
head_metadata_file=${main_metadata_dir}head_metadata.csv # testing
sample_list=test_sample_list.txt
temp_json_files_list=files.tmp
# find_comp_mutations.py files
PCM_data_file=${results_dir}${drug_of_interest}_PCM_data.txt
# filter_PCM.R files
PCM_model_results_file=${results_dir}${drug_of_interest}_PCM_model_results.csv
# clean_PCM files
PCM_merged_file=${results_dir}${drug_of_interest}_PCM_merged.csv
# comp_mut2res_mut.py files
PRM_stats_file=${results_dir}${drug_of_interest}_PRM_stats.csv
PRM_samples_file=${results_dir}${drug_of_interest}_PRM_samps.csv
samples_for_vcf_file=${results_dir}${drug_of_interest}_RM_samps.txt
binary_table_file=${results_dir}${drug_of_interest}_binary_table.csv
summary_file=${results_dir}${drug_of_interest}_summary.csv
# tree_pipeline.sh files
gvcf_file_suffix=.g.vcf.gz 
multi_samp_vcf=${vcf_dir}${drug_of_interest}.val.gt.g.vcf.gz
filt_multi_samp_vcf_file=${vcf_dir}${drug_of_interest}.filt.val.gt.g.vcf.gz
fasta_file=${fasta_dir}/${drug_of_interest}.filt.val.gt.g.snps.fa
newick_file=${newick_dir}${drug_of_interest}.filt.val.gt.g.snps.fa.treefile 
# comp_mut_tree.R
tree_png=${results_dir}${drug_of_interest}_tree.png

# ---------------------------------------------------
# Find all novel comp mutations for drug of interest
# ---------------------------------------------------

if [ ! -f ${PCM_data_file} ]; then
echo " --- PULLING COMPENSATORY MUTATIONS DATA FOR ${drug_of_interest} RUNNING python_scripts/find_PCM.py --- "
python python_scripts/find_PCM.py \
--metadata-file ${main_tb_metadata_file} \
--KCM-file ${KCM_file} \
--drug-of-interest ${drug_of_interest} \
--id-key ${id_key} \
--tbp-results ${tbp_results_dir} \
--outfile ${PCM_data_file}
fi

# ------------------------------------------------------------------------------------------------------------------------------------
# Run filter_PCM.R - run GLM models to see if DST is significantly predicted against presence of each mutation and lineage
# ------------------------------------------------------------------------------------------------------------------------------------

if [ ! -f ${PCM_model_results_file} ]; then
echo " --- FILTERING POTENTIAL NEW COMPENSATORY MUTATIONS - RUNNING r_scripts/filter_PCM.R --- "
Rscript r_scripts/filter_PCM.R \
--metadata_file ${main_tb_metadata_file} \
--PCM_data_file ${PCM_data_file} \
--drug_of_interest ${drug_of_interest} \
--outfile ${PCM_model_results_file} 
fi

# --------------------------------------------------------------------------------------------------------
# Concatenate output of filter_PCM.R with previous analyses of potential comp. mutations by TC
# --------------------------------------------------------------------------------------------------------

echo " --- CLEANING AND MERGING ${tc_file} AND ${PCM_model_results_file} - RUNNING clean_PCM.R"
Rscript r_scripts/clean_PCM.R \
--drug_of_interest ${drug_of_interest} \
--tc_file ${tc_file} \
--gn_results_file ${PCM_model_results_file} \
--KCM_file ${KCM_file} \
--outfile ${PCM_merged_file}
echo "--- output from clean_PCM.R - PCM for ${drug_of_interest} ---"
cat ${PCM_merged_file}
echo ""


# ----------------------------------------------------------------------------------------
# Filter potential novel compensatory mutations with tbprofiler critera 
# (n lineages, proportion of samples drug resistant, proportion of samples DST resistant)
# ----------------------------------------------------------------------------------------

echo " --- GETTING POTENTIAL NEW RESISTANCE MUTATIONS FOR ${drug_of_interest}; RUNNING python_scripts/comp_mut2res_mut.py --- "
python python_scripts/comp_mut2res_mut.py \
--drug-of-interest ${drug_of_interest} \
--PCM-file ${PCM_merged_file} \
--metadata-file ${main_tb_metadata_file} \
--tbdb-file ${tbdb_file} \
--drtypes-file ${dr_types_file} \
--KCM-file ${KCM_file} \
--tbprofiler-results-dir ${tbp_results_dir} \
--vars-exclude-file ${vars_exclude_file} \
--PRM-stats-file ${PRM_stats_file} \
--PRM-samples-file ${PRM_samples_file} \
--binary-table-file ${binary_table_file} \
--summary-file ${summary_file}

# ----------------------------
# Put the samples into a file
# ----------------------------

cat ${PRM_samples_file} | csvtk grep -f drug -p ${drug_of_interest} | csvtk cut -f wgs_id | tail -n+2 > ${samples_for_vcf_file}

# ------
# Trees
# ------

# Run tree pipeline to get newick tree file
if [ ! -f ${newick_file} ]; then
tree_pipeline.sh ${drug_of_interest} ${vcf_db} ${samples_for_vcf_file} ${vcf_dir} ${fasta_dir} ${newick_dir}
fi

# Plot tree and save as png
Rscript r_scripts/comp_mut_tree.R \
--tree_file ${newick_file} \
--metadata_file ${PRM_samples_file} \
--project_code ${drug_of_interest} \
--column 'drug' \
--outfile ${tree_png}

