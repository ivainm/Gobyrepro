#!/bin/bash
#SBATCH --partition=CPUQ  #use partition "short" for test runs
#SBATCH --account=nv-bio
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name="poplarge"
#SBATCH --mail-user=ivain.martinossi@ntnu.no
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err


cd ${SLURM_SUBMIT_DIR}
cd /cluster/home/ivainm/2.without_adapter/final_Aligned3

mkdir outpop_low_stringency

cd ..
export PATH=$PATH:/cluster/home/ivainm/programs/stacks/usr/bin   

# -p minimum number of populations a locus must be present in to be kept . 1= keep all
# -r minium percentage of individuals per pop to have the locus to keep it.
# -R minimum percentage accross pops. These should all disappear when I filter for missing data later anyways... I could be more stringent

#normal pipeline
populations -t 20 --structure --vcf --fstats -P ./final_Aligned3/catalogg_nomerge -M ../popmaps/pop_map_SORT_final.txt -O ./final_Aligned3/outpop_low_stringency -p 1 -r 0.1 -R 0.1

#for generating whitelist
#populations -t 20 --fasta-samples --vcf -P ./final_Aligned3/catalogg_nomerge_migrate -M ../popmaps/pop_map_SORT_migrate.txt -O ./final_Aligned3/outpop_migrate3 -p 6 -r 1. -R 1.
#for migrate input file if using whitelist
#populations -t 20 --fasta-samples --vcf -P ./final_Aligned3/catalogg_nomerge_migrate -M ../popmaps/pop_map_SORT_migrate.txt -O ./final_Aligned3/outpop_migrate4 -W ./final_Aligned3/outpop_migrate/whitelist.tsv -p 6 -r 1. -R 1.
#not using whitelist
#populations -t 20 --fasta-samples --vcf -P ./final_Aligned3/catalogg_nomerge_migrate -M ../popmaps/pop_map_SORT_migrate.txt -O ./final_Aligned3/outpop_migrate4 -p 3 -r 0.5 -R 0.5

#For identifying sex linked SNPs and then getting the consensus sequence just for these.
#mkdir outpop_sex
#Run normal pipeline first, then use the vcf to get the SNP ID (plink association test)
#Run population again from the same catalogue with whitelist of SNPs. The whitelist only needs to contain the vatalog locus ID. In the VCF file it should be the first element of the XX:XX:XX column
#populations -t 20 --fasta-samples --fasta-loci --vcf -P ./final_Aligned3/catalogg_nomerge -M ../popmaps/pop_map_SORT_final.txt -O ./final_Aligned3/outpop_sex -W ./final_Aligned3/outpop_sex/sex_candidates_whitelist.tsv -p 1 -r 0.1 -R 0.1
