#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 1 #number of compute units
#SBATCH -J post_trim_qc
#SBATCH -t 1:00:00 #max time

module load bioinfo-tools FastQC/0.11.9

#copy data to node local storage
SRC_DIR=/home/phennig/genome_analysis/qc_and_trim/rna_trim
OUT_DIR=/home/phennig/genome_analysis/qc_and_trim/trimmed_rna_fastqc_out/

mkdir $SNIC_TMP/fastqc_out

cp $SRC_DIR/*P.gz $SNIC_TMP/

#run code
fastqc $(ls $SNIC_TMP/*P.gz) -o $SNIC_TMP/fastqc_out

#copy results back
cp $SNIC_TMP/fastqc_out/* $OUT_DIR/
