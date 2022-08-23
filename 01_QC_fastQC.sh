#!/bin/bash -l

#SBATCH -J fastQC
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A xxx
#SBATCH -M snowy

#load modules:
module load bioinfo-tools
module load FastQC/0.11.9
module load MultiQC/1.9


input=$1 ##bam file
output_dir=$2

fastqc $input -o $output_dir 


