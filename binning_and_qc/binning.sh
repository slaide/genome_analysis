#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 4 #number of compute units
#SBATCH -J binning_qc
#SBATCH -t 4:00:00 #max time

module load bioinfo-tools MetaBat/2.12.1

metabat -t 1 -v -i ../assembly_and_qc/dna_assembly_combined/out_as/combined.contigs.fa -o binned_contigs
