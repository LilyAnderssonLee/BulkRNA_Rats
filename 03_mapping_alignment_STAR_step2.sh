#!/bin/bash -l
#SBATCH -J STAR-alignment
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 2-00:00:00
#SBATCH -A xxx
#SBATCH -M snowy

ulimit -c unlimited
set -eo pipefail

# load modules
module load bioinfo-tools
module load star/2.7.2b
ml samtools/1.9


#path to indexed genome
indexed_ref=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/star_mapping/star_index/
#path to trimmed fq
trimed_fq=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/quality_control/trim_galore
#fq.gz file name
read_fq=$1



star \
--runMode alignReads \
--genomeDir $indexed_ref \
--readFilesIn ${trimed_fq}/${read_fq}_trimmed.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix ${read_fq}_ \
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate


samtools index ${read_fq}_Aligned.sortedByCoord.out.bam
samtools flagstat ${read_fq}_Aligned.sortedByCoord.out.bam > ${read_fq}_STAR_flagstat.txt


