#!/bin/bash

# Run function findPhysLink with GEMMA v0.96 only for association analysis
# Copyright 2019 Yu Wan (wanyuac@gmail.com)
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 21 April 2019

#SBATCH --account=[your account name]
#SBATCH --job-name=GeneMatesLMM
#SBATCH --partition=[partition name]
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8192
#SBATCH --time=12:00:00

cd [workding directory]

module load [R model name]  # an optional step

export LANG="en_US.UTF-8"
export LC_COLLATE="en_US.UTF-8"
export LC_CTYPE="en_US.UTF-8"
export LC_MESSAGES="en_US.UTF-8"
export LC_MONETARY="en_US.UTF-8"
export LC_NUMERIC="en_US.UTF-8"
export LC_TIME="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"

Rscript --vanilla lmm.R
