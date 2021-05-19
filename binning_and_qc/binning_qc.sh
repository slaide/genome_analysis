#!/bin/bash -l
#SBATCH -A g2021012 -M snowy
#SBATCH -p core #compute units
#SBATCH -n 4 #number of compute units
#SBATCH -J binning_qc
#SBATCH -t 4:00:00 #max time

module load bioinfo-tools
module load CheckM/1.0.12

checkm lineage_wf -x fa --reduced_tree -t 4 -f checkm.out ./binning_out binning_qc_out
