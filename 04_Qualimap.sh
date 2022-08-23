#!/bin/bash -l

#SBATCH -J Qualimap
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -A xxx
#SBATCH -M snowy

ulimit -c unlimited
set -eo pipefail

ml bioinfo-tools
ml QualiMap/2.2.1

#By default, Qualimap will try to open a GUI to run Qualimap, so we need to run the unset DISPLAY command
unset DISPLAY

star_sorted_bam=$1
gtf=/crex/proj/xxx/private/Lili/analysis/bulkRNA/ratGenome_rn6/genes/rn6.ncbiRefSeq.gtf

qualimap rnaseq \
-outdir /crex1/proj/xxx/private/Lili/analysis/bulkRNA/qualimap/ \
-outfile $star_sorted_bam \
-outformat PDF \
-a proportional \
-bam /crex1/proj/xxx/private/Lili/analysis/bulkRNA/star_mapping/alignment/${star_sorted_bam}_Aligned.sortedByCoord.out.bam \
-gtf $gtf \
--java-mem-size=6G


