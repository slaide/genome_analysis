#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 1 #number of compute units
#SBATCH -J determine_expression_levels
#SBATCH -t 6:00:00 #max time

module load bioinfo-tools htseq/0.12.4
module load bioinfo-tools samtools/1.12

i=0

#only for those bins that fulfill the quality requirements from the binning qc post-processing (manually performed)
for bin_num in 18 10 21  3  7 16  5  9 24 13 20 14 22 17  4 15  2 ; do

echo "-------------------------------------------------------------------------------------"
echo "running for bin ${bin_num} (${i}/17)"
i=$((i+1))
date
echo "-------------------------------------------------------------------------------------"

# if output is already present, skip this bin
if ! test -s bin_${bin_num}_expression.csv ; then

# if input file for counting does not exist, create intermediate files and delete them afterwards
bamfilesconcatsorted="bin${bin_num}_concat_sorted.bam"
if ! test -s $bamfilesconcatsorted ; then

echo "concatenated and sorted mapped rna reads not present. creating."

echo "converting sam to bam files" ; date

sample1samfile="../rna_mapping/map_sample1_bin${bin_num}.sam"
sample2samfile="../rna_mapping/map_sample2_bin${bin_num}.sam"

# convert sam to bam files to save storage
sample1bamfile="sample1_bin${bin_num}.bam"
sample2bamfile="sample2_bin${bin_num}.bam"

nextcmd="samtools view -b -o $sample1bamfile $sample1samfile"
echo "running ${nextcmd}"
$nextcmd
nextcmd="samtools view -b -o $sample2bamfile $sample2samfile"
echo "running ${nextcmd}"
$nextcmd

echo "concatenating bam files" ; date
#concatenate both bam files into a single bam file
bamfilesconcat="bin${bin_num}_concat.bam"

nextcmd="samtools cat -o $bamfilesconcat $sample1bamfile $sample2bamfile"
echo "running ${nextcmd}"
$nextcmd

rm $sample1bamfile $sample2bamfile

echo "sorting concatenated bam file" ; date
# sorted concatenated bam file
#bamfilesconcatsorted="bin${bin_num}_concat_sorted.bam"
#if ! test -s $bamfilesconcatsorted ; then
samtools sort -o $bamfilesconcatsorted -O BAM $bamfilesconcat
rm $bamfilesconcat

nextcmd="samtools index ${bamfilesconcatsorted}"
echo "running ${nextcmd}"
$nextcmd

#end if sorted concat file did not exist
fi

# do the counting for concatenated sorted bam and gff file
gfffile="../annotation/prokka_out/bin_${bin_num}/bin_${bin_num}_nofa.gff"

# if ! test -s bin_${bin_num}_expression.csv ; then
echo "counting reads" ; date
# count paired end reads, sorted by position, write output to csv file, use CDS as feature type and use attribute with name ID as id
nextcmd="htseq-count -d , -c bin_${bin_num}_expression.csv -i ID -t CDS -r pos $bamfilesconcatsorted $gfffile"
echo "running ${nextcmd}"
$nextcmd

# if output already present
else
echo "expression analysis output already present. skipping."

# end if output already present
fi

done
