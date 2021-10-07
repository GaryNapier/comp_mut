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
tbp_results_dir=~jody/tbprofiler_tests/results/

# -----------------------------------------------------------
# Run tbprofiler_template.py to find novel mutations in katG
# -----------------------------------------------------------


# python python_scripts/tbprofiler_template_testing.py --sample_list_file test_sample_list --dir ../pakistan/tbprofiler_pakistan_results/json

python python_scripts/tbprofiler_template.py --dir ${tbp_results_dir}

