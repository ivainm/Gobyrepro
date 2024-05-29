##############################################################
README for the Gobyrepro project

#############################################
The project involves two distinct analysis:
1. Analysis involving ddRAD sequencing data. 
	-Read processing and filtering, quality checks
	-Alignment to reference genome
	-Using the Stacks pipeline to find SNPs and produce a VCF file
	-Using VCFtools and Plink to filter the SNPs 
	-Use Plink and R to do some bioinformatics analysis (mainly population structure)
	
2. Analysis involving field collected data and SNP data from analysis 1 to produce manuscript figures and perform statsitical tests:
	-All done in R
	-Different types of data: phenotypes, population census, environmental data, SNP data.
	
	
	
	
############################################
Part 1
Slurm code were run on the IDUN cluster computer of NTNU. Bash and R scripts on local.

Data files             #description
-------------------------------------
###repository address to raw data
pop.conv.tsv           #file assigning individuals to their population, in the format required for the radiator package (Fst calculation in R)
R90_maf.bed/bim/fam    #3 files constituting the bed file obtained from the filtered vcf file, in order to calcualte Fst with the radiator package in R


Code files             #description									
-------------------------------------
processrad.sh          #demultiplexing the raw reads and filtering for adapter content. Removing suspicious technical replicates an merging remaining ones.
align2.slurm           #align demultiplexed ddRAD reads to the reference genome of P.minutus. Index aligned reads. Run som quality checks with samtools.
gstacks.slurm          #firt step of the Stacks pipeline - SNP catalogue creation from aligned reads.
populations.slurm      #second step of the Stacks pipeline - VCF file output.
final_vcf.sh           #from the vcf file, filter using vcftools and plink, run pca, export a bedfile for use in R, calculate individual level heterozygosity.
R_popgen.R             #Different step of the analysis done in R. Plots to identify outliers, Fst calculations with radiator package.







###########################################
Part 2

Data files             #description
-------------------------------------
ALLmod.csv             #Sea surface temperature data from Nordkyst800 model for the studied populations
nestdata.csv           #Field data collected in artificial nests (egg features, phenotypes of guarding males)
OSRdata.csv            #Field data of population census from observed individuals
phenodata.csv          #Field data of fish phenotypes of captured individuals
sexrpheno.csv          #Census data from captured individuals

Allfit.txt             #Fitted sea surface temperature curves for the studied populations
Coord.txt              #Coordinates of the study populations
devtemp.txt            #Experimental data of egg development speed at different water temperatures, from previous study
Dist.txt               #Coastal distance in km between the different populations
env_info.txt           #Modelled developmental time available for each population during one season
out_all.txt            #List of outlier individuals to be removed when analysing genomic data
pop_map_final.txt      #Information about individuals sequenced
primaryprod.txt        #Estimated beginning and end dates of primary productivity in the studied populations 


Code files             #description
-------------------------------------
final_script_REPRO.R   #Using the data files above, produces all the statistical tests and plots of the manuscript.
