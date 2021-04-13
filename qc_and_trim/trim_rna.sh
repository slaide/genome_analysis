#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 4 #number of compute units
#SBATCH -J trim_rna
#SBATCH -t 2:00:00 #max time

module load bioinfo-tools trimmomatic/0.36

#copy data to node local storage
SRC_DIR=/home/phennig/genome_analysis/raw_data
OUT_DIR=/home/phennig/genome_analysis/qc_and_trim

cp $SRC_DIR/RNA/* $SNIC_TMP

#run code

trimmomatic PE -threads 4 -basein $SNIC_TMP/SRR4342137.1.fastq.gz -baseout $SNIC_TMP/SRR4342137_filtered.fastaq.gz ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa:2:30:10 \
                LEADING:6 TRAILING:6 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:30
trimmomatic PE -threads 4 -basein $SNIC_TMP/SRR4342139.1.fastq.gz -baseout $SNIC_TMP/SRR4342139_filtered.fastaq.gz ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa:2:30:10 \
                LEADING:6 TRAILING:6 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:30

#copy results back
cp $SNIC_TMP/*filtered* $OUT_DIR/
