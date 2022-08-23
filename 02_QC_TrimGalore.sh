#!/bin/bash -l

#SBATCH -J Trim
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -A xxx
#SBATCH -M snowy

module load bioinfo-tools
module load FastQC/0.11.9
module load MultiQC/1.9
#module load spider
module load TrimGalore/0.6.1

#load fastq.gz file
fq=$1
output_dir=$2

trim_galore $fq --quality 20 --illumina --phred33 --length 20 --gzip --output_dir $output_dir --fastqc 

