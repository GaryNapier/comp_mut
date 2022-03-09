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
newick_dir=newick/
vcf_db=~/vcf/
vcf_dir=vcf/

# ------
# Files
# ------

# Existing files
main_tb_metadata_file=${main_metadata_dir}tb_data_18_02_2021.csv
known_comp_mut_file=${database_dir}compensatory_mutations.csv 
tbdb_file=${tbdb_dir}tbdb.csv
dr_types_file=${database_dir}dr_types.json
vars_exclude_file=${local_metadata_dir}var_exclude_comp_mut.csv

# Created files
head_metadata_file=${main_metadata_dir}head_metadata.csv # testing
sample_list=test_sample_list.txt
temp_json_files_list=files.tmp
# find_comp_mutations.py files
novel_comp_mut_data_file=${results_dir}${drug_of_interest}_novel_comp_mut_data.txt
# filter_comp_mut.R files
novel_comp_mut_model_results_file=${results_dir}${drug_of_interest}_novel_comp_mut_model_results.csv
# comp_mut2res_mut.py files
potential_res_mut_stats_file=${results_dir}potential_res_mut_stats.csv
potential_res_mut_samples_file=${results_dir}potential_res_mut_samps.csv
samples_for_vcf_file=${results_dir}${drug_of_interest}_res_mut_samps.txt

# variant_calling_and_concat_gvcfs.sh
gvcf_file_suffix=.g.vcf.gz 
# variant_filtering.sh
multi_samp_vcf=${vcf_dir}${drug_of_interest}.val.gt.g.vcf.gz
filt_multi_samp_vcf_file=${vcf_dir}${drug_of_interest}.filt.val.gt.g.vcf.gz
# vcf2fasta.sh
# iqtree.sh
fasta_file=${fasta_dir}/${drug_of_interest}.filt.val.gt.g.snps.fa


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

# ----------------------------------------------------------------------------------------
# Filter potential novel compensatory mutations with tbprofiler critera 
# (n lineages, proportion of samples drug resistant, proportion of samples DST resistant)
# ----------------------------------------------------------------------------------------

echo " --- GETTING POTENTIAL NEW RESISTANCE MUTATIONS FOR ${drug_of_interest}; RUNNING python_scripts/comp_mut2res_mut.py --- "
python python_scripts/comp_mut2res_mut.py \
--drug-of-interest ${drug_of_interest} \
--potential-comp-mut-file ${novel_comp_mut_model_results_file} \
--metadata-file ${main_tb_metadata_file} \
--tbdb-file ${tbdb_file} \
--drtypes-file ${dr_types_file} \
--known-comp-mut-file ${known_comp_mut_file} \
--tbprofiler-results-dir ${tbp_results_dir} \
--vars-exclude-file ${vars_exclude_file} \
--potential-res-mut-stats-file ${potential_res_mut_stats_file} \
--potential-res-mut-samples-file ${potential_res_mut_samples_file}

# ----------------------------
# Put the samples into a file
# ----------------------------

cat ${potential_res_mut_samples_file} | csvtk grep -f drug -p isoniazid | csvtk cut -f samples | tail -n+2 > ${samples_for_vcf_file}

# ------
# Trees
# ------

tree_pipeline.sh ${drug_of_interest} ${vcf_db} ${samples_for_vcf_file} ${vcf_dir} ${fasta_dir} ${newick_dir}
