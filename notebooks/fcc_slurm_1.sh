#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lir9@ccf.org
#SBATCH --job-name=vqec_p_0_2_u_0_1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH -p defq
#SBATCH --mem=100000
#SBATCH -o fcc_vqec_p_0_2_u_0_1_%A.out
#SBATCH -e fcc_vqec_p_0_2_u_0_1_%A.err
python fcc_vqec_sim_1.py
