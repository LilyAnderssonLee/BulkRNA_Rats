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
count_matrix=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/rawCounts/featureCounts/rats_C9N2AANXX_7_20160809B_20160809_featureCounts.txt
star_bam=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/star_mapping/alignment/*Aligned.sortedByCoord.out.bam


featureCounts -T 4 \
-a $anno \
-o $count_matrix \
-g gene_id \
-O \
-t exon \
--fraction \
$star_bam

 

# Strandness

#  -s <int or string>  Perform strand-specific read counting. A single integer
#                      value (applied to all input files) or a string of comma-
#                      separated values (applied to each corresponding input
#                      file) should be provided. Possible values include:
#                      0 (unstranded), 1 (stranded) and 2 (reversely stranded).
#                      Default value is 0 (ie. unstranded read counting carried
#                      out for all input files).


#When in doubt, run all 3 on one sample. If the counts from yes and reverse are roughly equal for most genes then the dataset is unstranded (i.e., no is the correct setting). 
#If either yes or reverse produces much higher counts than the other then the appropriate setting is the one giving the higher counts. This will pretty much always be reverse these days. 
#can also just ask the people who did the library prep. what kit they used or whether it was dUTP-based.