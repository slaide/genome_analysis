#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 2 #number of compute units
#SBATCH -J metagenome_assembly_and_qc_w_combined_data
#SBATCH -t 12:00:00 #max time

module load bioinfo-tools megahit/1.2.9

SRC_DIR=/home/phennig/genome_analysis/raw_data/DNA
OUT_DIR=/home/phennig/genome_analysis/assembly_and_qc/dna_assembly_combined

cp -r $SRC_DIR/ $SNIC_TMP/in/

#run code
echo "starting assembly"
date

megahit --k-min 65 --k-max 105 --k-step 10 -t 2 --kmin-1pass \
-1 $SNIC_TMP/in/SRR4342129_1.paired.trimmed.fastq.gz,$SNIC_TMP/in/SRR4342133_1.paired.trimmed.fastq.gz \
-2 $SNIC_TMP/in/SRR4342129_2.paired.trimmed.fastq.gz,$SNIC_TMP/in/SRR4342133_2.paired.trimmed.fastq.gz  \
-o $SNIC_TMP/out_as --out-prefix combined

echo "done with assembly"
date

#copy results back
echo "copy assembly data out"
cp -r $SNIC_TMP/out_as $OUT_DIR/

echo "starting qc"
date

module load bioinfo-tools quast/5.0.2

SRC_DIR=/home/phennig/genome_analysis/assembly_and_qc/dna_assembly_combined/out_as
OUT_DIR=/home/phennig/genome_analysis/assembly_and_qc/dna_qc_combined

mkdir $SNIC_TMP/out_qc

cp -r $SRC_DIR/ $SNIC_TMP/in/

#run code
echo "eval assembly"
date

quast.py -t 2 -o $SNIC_TMP/out_qc $SNIC_TMP/in/combined.contigs.fa

echo "eval done"
date

#copy results back
cp -r $SNIC_TMP/out_qc $OUT_DIR/

echo "done"
