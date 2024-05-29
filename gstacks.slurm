#!/bin/bash
#SBATCH --partition=CPUQ  #use partition "short" for test runs
#SBATCH --account=nv-bio
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name="gstacks2"
#SBATCH --mail-user=ivain.martinossi@ntnu.no
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err


cd ${SLURM_SUBMIT_DIR}
cd /cluster/home/ivainm/2.without_adapter/final_Aligned3

mkdir catalogg_nomerge
#mkdir catalogg_nomerge_migrate
cd ..

export PATH=$PATH:/cluster/home/ivainm/programs/stacks/usr/bin    

gstacks -t 20 -I ./final_Aligned3/sorted -M ../popmaps/pop_map_SORT_final.txt -O ./final_Aligned3/catalogg_nomerge  #bamfiles have to be sorted or error message
#gstacks -t 20 -I ./final_Aligned3/sorted -M ../popmaps/pop_map_SORT_migrate.txt -O ./final_Aligned3/catalogg_nomerge_migrate 