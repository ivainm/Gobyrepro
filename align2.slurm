#!/bin/bash
#SBATCH --partition=CPUQ  #use partition "short" for test runs
#SBATCH --account=nv-bio
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="align4"
#SBATCH --mail-user=ivain.martinossi@ntnu.no
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

cd ${SLURM_SUBMIT_DIR}

#This slurm code does the following tasks
#Align demultiplexed reads to the P.minutus reference genome
#convert the sam files to bam files
#sort the bam files
#index the bam files
#produce stats report for multiqc to assess alignment success
#produce a multiqc report


echo "we are running from this directory: $SLURM_SUBMIT_DIR"
echo " the name of the job is: $SLURM_JOB_NAME"
echo "Th job ID is $SLURM_JOB_ID"
echo "The job was run on these nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "We are using $SLURM_CPUS_ON_NODE cores"
echo "We are using $SLURM_CPUS_ON_NODE cores per node"
echo "Total of $SLURM_NTASKS cores"

cd /cluster/home/ivainm
cd 2.without_adapter
mkdir final_Aligned3
cd final_Aligned3
mkdir sorted
cd sorted
mkdir stats
cd ..
cd ..

export PATH=$PATH:/cluster/home/ivainm/programs/samtools/usr/bin/ #only when running on IDUN

for rep in 1 2 3
do
for var in $(seq -f "%02g" 1 180)
do
~/programs/bwamem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 20 ~/indexbwa/bwaindex ./final_merged2/DYN-$var".1".fq.gz ./final_merged2/DYN-$var".2".fq.gz | samtools view --threads 20 -S -b -h - > ./final_Aligned3/ALI-$var.bam
samtools sort ./final_Aligned3/ALI-$var.bam -o ./final_Aligned3/sorted/SORT-$var.bam    #sort in order of ref sequence
samtools index ./final_Aligned3/sorted/SORT-$var.bam	
samtools stats --threads 20 ./final_Aligned3/sorted/SORT-$var.bam > ./final_Aligned3/sorted/stats/SORT-$var.stats #create stats file for multiQC
done
done

module load Python/3.10.4-GCCcore-11.3.0
pip install multiqc  
multiqc ./final_Aligned3/sorted/stats
module purge
