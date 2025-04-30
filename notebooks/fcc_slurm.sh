#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lir9@ccf.org
#SBATCH --job-name=fcc_vqec
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -p defq
#SBATCH --mem=100000
#SBATCH -o fcc_vqec_%A.out
#SBATCH -e fcc_vqec_%A.err
python fcc_vqec_sim.py
