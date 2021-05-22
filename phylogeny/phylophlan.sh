#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 2 #number of compute units
#SBATCH -J phylophlan
#SBATCH -t 10:00:00 #max time

module load bioinfo-tools phylophlan muscle FastTree usearch/5.2.32

# copy prokka (binning) output files to local directory
# but not all bins! these listed below are those that were manully selected (those not present were omitted due to unfit metrics)
#mkdir -p input/prokka_bins_faa
#for bin_num in 18 10 21  3  7 16  5  9 24 13 20 14 22 17  4 15  2 ; do
#target="../annotation/prokka_out/bin_${bin_num}/bin_${bin_num}.faa"
#link="input/prokka_bins_faa/bin_${bin_num}.faa"
#cp $target $link
#done

phylophlan.py --nproc 2 -t -i prokka_bins_faa
echo "___________________________________________________________"
echo "done"
