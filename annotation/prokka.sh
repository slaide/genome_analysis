#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 2 #number of compute units
#SBATCH -J annotation
#SBATCH -t 4:00:00 #max time

module load bioinfo-tools prokka/1.45-5b58020

loop_num=1

#manually specified bin nums with contamination<10% and completeness>10% (17/26 bins)
for bin_num in 18 10 21  3  7 16  5  9 24 13 20 14 22 17  4 15  2
do

echo "__________________________________________________________________________"
echo "running loop ${loop_num}/17 for bin ${bin_num}"
echo "__________________________________________________________________________"
prokka --addgenes --metagenome --outdir ./prokka_out/bin_${bin_num} --cpus 2 --prefix bin_${bin_num} ../binning_and_qc/binning_out/binned_contigs_${bin_num}.fa 

loop_num=$((loop_num+1))

done
