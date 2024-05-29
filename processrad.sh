#!/bin/bash

#Gobyrepro project
#This code is for demultiplexing raw ddRAD reads and filtering adapter content
#First chunk does not merge technical replicates (for the replicate filtering stage)
#Second step merges replicates after removing suspicious technical replicates (bad apples)

#Chunk 1
cd ~/realdir
mkdir newdemult_nomerge
ulimit -n 3500
#Demultiplex and remove expected adapter content
process_radtags -r -P -q -p ./rawreads/ddRAD-DYN01/ -o ./newdemult_nomerge/ -b ./rawreads/barcodes_DYN01.txt \
--inline_inline --renz_1 pstI --renz_2 mseI --threads 12 --adapter-mm 3 --adapter-1 ACACTCTTTCCCTACACGACGCTCTTCCGATCT --adapter-2 AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCAGAACAA
 
process_radtags -r -P -q -p ./rawreads/ddRAD-DYN02/ -o ./newdemult_nomerge/ -b ./rawreads/barcodes_DYN02.txt \
 --inline_inline --renz_1 pstI --renz_2 mseI --threads 12 --adapter-mm 3 --adapter-1 ACACTCTTTCCCTACACGACGCTCTTCCGATCT --adapter-2 AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCAGAACAA
 
mkdir nomerge_cutadapt
 #remove adapter content with generic illumina sequence
for rep in 1 2 3
do for var in $(seq -f "%02g" 1 180)
do cutadapt -j 12 -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
 --discard-trimmed \
 -o ./nomerge_cutadapt/DYN-$var"_0"$rep".1".fq.gz -p ./nomerge_cutadapt/DYN-$var"_0"$rep".2".fq.gz ./newdemult_nomerge/DYN-$var"_0"$rep".1".fq.gz ./newdemult_nomerge/DYN-$var"_0"$rep".2".fq.gz 
done
done

cd nomerge_cutadapt
#get some basic quality checks from fastqc
mkdir stats
fastqc -t 12 --outdir stats ./DYN-*


#Chunk 2
#Merge the replicates after removing the bad apples
cd demult_final
mkdir final_merged2

#I need to include an if statement to check if replicates belong to the bad apples before merging
#just remember that the file containing the names to exclude must be converted to unix before use!

cd ../nomerge_cutadapt

for pair in 1 2
do for var in $(seq -f "%02g" 1 180)
do printf 'DYN-'$var'_01.'$pair'.fq.gz\nDYN-'$var'_02.'$pair'.fq.gz\nDYN-'$var'_03.'$pair'.fq.gz' | grep -Fvf ../demult_final/remove_reps.txt | xargs cat > ../demult_final/final_merged2/DYN-$var.$pair.fq.gz
done
done


#for pair in 1 2
#do for var in $(seq -f "%02g" 1 180)
#do cat ../nomerge_cutadapt/DYN-$var"_01".$pair.fq.gz ../nomerge_cutadapt/DYN-$var"_02".$pair.fq.gz ../nomerge_cutadapt/DYN-$var"_03".$pair.fq.gz > final_merged/DYN-$var.$pair.fq.gz
#done
#done

#for pair in 1 2
#do for var in $(seq -f "%02g" 96 180)
#do cat newdemult2_adp_mm2/DYN-$var"_01".$pair.fq.gz newdemult2_adp_mm2/DYN-$var"_02".$pair.fq.gz newdemult2_adp_mm2/DYN-$var"_03".$pair.fq.gz > raw_merged_adp_mm2/DYN-$var.$pair.fq.gz
#done
#done
