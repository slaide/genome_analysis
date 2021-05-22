#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 2 #number of compute units
#SBATCH -J rna_mapping_bwa_mem
#SBATCH -t 2:00:00 #max time

module load bioinfo-tools bwa/0.7.17 samtools/1.12

i=1

for bin_num in 18 10 21  3  7 16  5  9 24 13 20 14 22 17  4 15  2 ; do

echo "starting rna mapping for bin ${bin_num} (${i}/17)"
date

bin_dna_file="../binning_and_qc/binning_out/binned_contigs_${bin_num}.fa"

bwa index ${bin_dna_file}

outfile1="map_sample1_bin${bin_num}.sam"
outfile2="map_sample2_bin${bin_num}.sam"

bwa mem -o $outfile1 -t 2 ${bin_dna_file} ../qc_and_trim/rna_trim/SRR4342137_filtered.fastaq_1P ../qc_and_trim/rna_trim/SRR4342137_filtered.fastaq_2P
bwa mem -o $outfile2 -t 2 ${bin_dna_file} ../qc_and_trim/rna_trim/SRR4342139_filtered.fastaq_1P ../qc_and_trim/rna_trim/SRR4342139_filtered.fastaq_2P

i=$((i+1))

done

echo "done"
