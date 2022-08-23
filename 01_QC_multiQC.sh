#!/bin/bash -l

#SBATCH -J multiQC
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A sxxx
#SBATCH -M snowy

#load modules:
module load bioinfo-tools
module load FastQC/0.11.9
module load MultiQC/1.9


input=$1 #the folder stores fastqc reports
output_dir=$2

multiqc $input -o $output_dir

