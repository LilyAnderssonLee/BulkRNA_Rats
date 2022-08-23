#!/bin/bash -l

#SBATCH -J htseq
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -A xxx
#SBATCH -M snowy

ml bioinfo-tools
ml htseq/0.9.1


star_bam_path=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/star_mapping/alignment
gtf=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/ratGenome_rn6/genes/rn6.ncbiRefSeq.gtf
bam_name=$1

htseq-count -s no -r pos -t exon -i gene_id -f bam ${star_bam_path}/${bam_name}_Aligned.sortedByCoord.out.bam $gtf >${bam_name}_htseq_counts


