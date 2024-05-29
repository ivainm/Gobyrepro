#!/bin/bash
#In this script I analyze SNPs from non merged replicates
#The main point is to identify if the outlier individuals have consistently weird replicates or if
#some replicate can be isolated and excluded and the rest kept so I don't have to remove the whole individuals
#I will also look more generally at replicate consistency across all individuals

cd ~/realdir/vcfs/2.without_adapter/nomerge

#Do I need to separate the replicates before filtration steps are applied?
#For example could it affect maf and frequency filters in a weird way? Frequencies and means should not be affected only absolute counts
#this script starts with the output of the "population" step of Stacks
vcftools --vcf populations.snps.vcf --min-alleles 2 --max-alleles 2 --recode --out maxbiallel
vcftools --vcf maxbiallel.recode.vcf --min-meanDP 5  --max-meanDP 30 --recode --out maxbiDP  
vcftools --vcf maxbiDP.recode.vcf --maf 0.05 --recode --out maxbiDPmaf5 #maf 0.05 instead of mac 3
vcftools --vcf maxbiDPmaf5.recode.vcf --max-missing 0.75 --recode --out maxbiDPmaf5_75 #90% missingness

#With the same stringency as for the merged data I get zero SNPs. Somehow there must be a lot of SNPs that are 
#only in one or two replicates so that the frequency filters are more stringent when it's not merged
#also, for the SNPs that are in all replicates, the depth filtering should also be different

vcftools  --vcf maxbiDPmaf5_75.recode.vcf --missing-indv --out missing_75 #check missingness per individual to see if some should be excluded

#I only keep individuals that have less than 75% missing sites on the previous file
chmod +rwx miss_out_75.txt
chmod +rwx miss_out_75_vcf.txt

#redo the depth, maf and missingness filtering, with those outliers removed
vcftools --vcf maxbiallel.recode.vcf --remove miss_out_75_vcf.txt  --max-meanDP 30 --recode --out maxbi_out  
vcftools --vcf maxbi_out.recode.vcf --min-meanDP 5  --max-meanDP 30 --recode --out maxbiDP_out  
vcftools --vcf maxbiDP_out.recode.vcf --maf 0.05 --recode --out maxbiDPmaf5_out #maf 0.05 instead of mac 3
vcftools --vcf maxbiDPmaf5_out.recode.vcf --max-missing 0.90 --recode --out maxbiDPmaf5_out90 #90% missingness


#bed file, HW filtering, LD pruning, pca, ped file for R, with outliers of missing data removed
~/plink_files/plink --vcf maxbiDPmaf5_out90.recode.vcf --remove miss_out_75.txt --make-bed --allow-extra-chr --out plink90
~/plink_files/plink --bfile plink90 --hwe 1e-10 midp include-nonctrl --nonfounders --make-bed --out plink90HW --allow-extra-chr
~/plink_files/plink --bfile plink90HW --indep-pairwise 50 5 0.5 --out plink90HWLD05 --allow-extra-chr
~/plink_files/plink --bfile plink90HW --extract plink90HWLD05.prune.in --pca --allow-extra-chr --out pca90_nomerge_out
~/plink_files/plink --bfile plink90HW --extract plink90HWLD05.prune.in --allow-extra-chr --recode --out plink90HW_R_out
