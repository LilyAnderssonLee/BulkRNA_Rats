#!/bin/bash -l

ml bioinfo-tools
ml MultiQC/1.9

log_file_dir=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/MultiQC_summary/log_files

multiqc --interactive $log_file_dir
