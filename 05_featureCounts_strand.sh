#!/bin/bash -l

#SBATCH -J featureCounts
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1-00:00:00
#SBATCH -A xxx
#SBATCH -M snowy

ulimit -c unlimited
set -eo pipefail

ml bioinfo-tools
ml subread/2.0.0


anno=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/ratGenome_rn6/genes/rn6.ncbiRefSeq.gtf
count_matrix=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/rawCounts/featureCounts/rats_C9N2AANXX_7_20160809B_20160809_featureCounts_Stranded.txt
star_bam=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/star_mapping/alignment/*Aligned.sortedByCoord.out.bam


featureCounts -T 4 \
-a $anno \
-o $count_matrix \
-g gene_id \
-O \
-t exon \
-s 1 \
--fraction \
$star_bam
