#jobid 4415399

#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 1 #number of compute units
#SBATCH -J pre_trim_qc
#SBATCH -t 1:00:00 #max time

module load bioinfo-tools FastQC/0.11.9

#copy data to node local storage
SRC_DIR=/home/phennig/genome_analysis/raw_data
OUT_DIR=/home/phennig/genome_analysis/qc_and_trim

mkdir $SNIC_TMP/DNA
mkdir $SNIC_TMP/RNA
cp $SRC_DIR/DNA/* $SNIC_TMP/DNA/
cp $SRC_DIR/RNA/* $SNIC_TMP/RNA/

#run code
mkdir $SNIC_TMP/dna_qc_fastqc_out
mkdir $SNIC_TMP/rna_qc_fastqc_out

fastqc $(ls $SNIC_TMP/DNA/*) -o $SNIC_TMP/dna_qc_fastqc_out
fastqc $(ls $SNIC_TMP/RNA/*) -o $SNIC_TMP/rna_qc_fastqc_out

#copy results back
cp -r $SNIC_TMP/dna_qc_fastqc_out $OUT_DIR/dna_qc_fastqc_out
cp -r $SNIC_TMP/rna_qc_fastqc_out $OUT_DIR/rna_qc_fastqc_out
