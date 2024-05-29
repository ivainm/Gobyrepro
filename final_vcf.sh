#!/bin/bash
#Gobyrepro project
#In this script I filter the vcf file from stacks, using vcftools and plink

cd ~/realdir/vcfs/2.without_adapter/final2

#First, I want to check what would be a reasonnable minimum depth to use for filtering. 
vcftools --vcf populations.snps.vcf --site-depth --out SNP_DP

#I start with the most basic filtration steps
vcftools --vcf populations.snps.vcf --min-alleles 2 --max-alleles 2 --recode --out maxbiallel
vcftools --vcf maxbiallel.recode.vcf --max-missing 0.90 --recode --out maxmiss90 #90% missingness

#vcftools --vcf maxmiss90.recode.vcf --min-meanDP 8  --max-meanDP 30 --recode --out maxmiss90_DP8_30  
#vcftools --vcf maxmiss90.recode.vcf --min-meanDP 10  --max-meanDP 30 --recode --out maxmiss90_DP10_30  
#vcftools --vcf maxmiss90.recode.vcf --min-meanDP 10  --max-meanDP 20 --recode --out maxmiss90_DP10_20 

vcftools --vcf maxmiss90_DP10_20.recode.vcf --maf 0.05 --recode --out DP1020maf5 #maf 0.05 instead of mac 3
#vcftools --vcf maxmiss90_DP10_20.recode.vcf --mac 2 --recode --out DP1020mac2  #mac 2 instead of maf 0.05
#the maf filtered file will be used to do neutral population structure, PCA, ADMIXTURE, Fst

vcftools  --vcf DP1020maf5.recode.vcf --missing-indv --out missingmaf5_90 #check missingness per individual to see if some should be excluded
#vcftools  --vcf DP1020mac2.recode.vcf --missing-indv --out missingmac3_90

#I follow with one step in R to visualize missingness and if some individuals should be excluded I export a file with their names.
#I only keep individuals that have less than 40% missing sites on the previous file
chmod +rwx miss_out_40.txt
chmod +rwx miss_out_40_vcf.txt
chmod +rwx out_all.txt

#Make pca in plink, output bedfile to use in R for admixture plot, and output heterozygosity
#to remove outliers: miss_out_40.txt is the list of individuals with >40% missingness. out_all.txt also includes the admixture and hz outliers.
#bed file, HW filtering, LD pruning, pca, ped file for admixture in R
~/plink_files/plink --vcf DP1020maf5.recode.vcf --remove out_all.txt --make-bed --allow-extra-chr --out plink90maf
~/plink_files/plink --bfile plink90maf --hwe 1e-10 midp include-nonctrl --nonfounders --make-bed --out plink90HWmaf --allow-extra-chr
~/plink_files/plink --bfile plink90HWmaf --indep-pairwise 50 5 0.5 --out plink90HWLD05maf --allow-extra-chr
#rewrite vcf version for later use e.g. dapc in R90_mac
~/plink_files/plink --bfile plink90HWmaf --extract plink90HWLD05maf.prune.in --recode vcf-iid --out plink90HWLD05mafR --allow-extra-chr
~/plink_files/plink --bfile plink90HWmaf --extract plink90HWLD05maf.prune.in --pca --allow-extra-chr --out pca90_maf  #pca
~/plink_files/plink --bfile plink90HWmaf --extract plink90HWLD05maf.prune.in --allow-extra-chr --recode --out plink90HW_R_maf #bed file export for use in R
~/plink_files/plink --bfile plink90HWmaf --extract plink90HWLD05maf.prune.in --het --allow-extra-chr --out het90_maf #heterozigosity

#~/plink_files/plink --bfile plink90HWmaf --extract plink90HWLD05maf.prune.in --make-bed --update-ids ../update_pops.txt --allow-extra-chr --out R90_maf
#~/plink_files/plink --bfile plink90HWmaf --extract plink90HWLD05maf.prune.in --make-bed --update-ids ../update_pops_season.txt --allow-extra-chr --out R90_maf_season

#Make different subsets for pcas
#Just Helligvaer and Ringstad, to check if 2yo breeders are different from 1yo breeders (we don't have unsexed) so far age is estimated from size but could be checked by Tilde if promissing
#~/plink_files/plink --bfile plink90HWmaf --within ../phenotypes.txt --mwithin 2 --keep-cluster-names Ringstad Helligvaer --extract plink90HWLD05maf.prune.in --pca --allow-extra-chr --out pca09maf_HellRig
#Just north pops, removing Kberg and Arendal
#~/plink_files/plink --bfile plink90HWmaf --within ../phenotypes.txt --mwithin 2 --keep-cluster-names Austevoll Hitra Ringstad Helligvaer --extract plink90HWLD05maf.prune.in --pca --allow-extra-chr --out pca09maf_North


# Association of allele frequencies with "age" as determined by size threshold in the two northern populations
#For such a test it is recommended not to LD prune the data and not to HWE filter either.
#~/plink_files/plink --bfile plink90maf --within ../phenoage.txt --mwithin 1 --keep-cluster-names Ringstad Helligvaer --pheno ../phenoage.txt --pheno-name AGE --assoc fisher --allow-extra-chr --allow-no-sex --out assocage

#Rerun the asociation test with sex now that the dataset has been filtered for outlier individuals etc
#keep filtration not too stringent in that case
#vcftools --vcf maxbiallel.recode.vcf --max-missing 0.75 --recode --out maxmiss75 #75% missingness
#vcftools --vcf maxmiss75.recode.vcf --min-meanDP 8  --max-meanDP 30 --recode --out maxmiss75_DP8_30
#vcftools --vcf maxmiss75_DP8_30.recode.vcf --mac 2 --recode --out DP75830mac2
#~/plink_files/plink --vcf DP75830mac2.recode.vcf --remove out_all.txt --make-bed --allow-extra-chr --out plink75mac
#let's go!
#~/plink_files/plink --bfile plink75mac --pheno ../phenotypes.txt --pheno-name SEX --assoc fisher --allow-extra-chr --allow-no-sex --out assocsex
#Visual insepction in R gives a p-vlue threshold of 10^-13 to select SNPs that look like good sex-linked candidates. I get 27 SNPs.
#In R I can export  list of the locus ID from the vcf that should allow me to whitelist them in Stacks/population




#Neutral markers for pca in R to give population levels PC loadings_ GIVES POPULATION LEVEL FREQUENCIES_ stringent HWE filtering and LD pruning
#~/plink_files/plink --bfile plink90maf --hwe 1e-5 midp include-nonctrl --nonfounders --make-bed --out plink90HWmaf_neutral --allow-extra-chr
#~/plink_files/plink --bfile plink90HWmaf_neutral --indep-pairwise 50 5 0.1 --out plink90HWLD01maf_neutral --allow-extra-chr
#~/plink_files/plink --bfile plink90HWmaf_neutral --extract plink90HWLD01maf_neutral.prune.in --within ./phenotypes.txt --mwithin 2 --freq --out plink09_pop_neutral --allow-no-sex --allow-extra-chr

#-----------------------------------------------
#Version with mac instead of maf filtering
#~/plink_files/plink --vcf DP1020mac2.recode.vcf --remove out_all.txt --make-bed --allow-extra-chr --out plink90mac
#~/plink_files/plink --bfile plink90mac --hwe 1e-10 midp include-nonctrl --nonfounders --make-bed --out plink90HWmac --allow-extra-chr
#~/plink_files/plink --bfile plink90HWmac --indep-pairwise 50 5 0.5 --out plink90HWLD05mac --allow-extra-chr
#~/plink_files/plink --bfile plink90HWmac --extract plink90HWLD05mac.prune.in --pca --allow-extra-chr --out pca90_mac
#~/plink_files/plink --bfile plink90HWmac --extract plink90HWLD05mac.prune.in --allow-extra-chr --recode --out plink90HW_R_mac

#~/plink_files/plink --bfile plink90HWmac --extract plink90HWLD05mac.prune.in --make-bed --update-ids ../update_pops.txt --allow-extra-chr --out R90_mac


