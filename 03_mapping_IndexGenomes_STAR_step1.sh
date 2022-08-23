#!/bin/bash -l
#SBATCH -J STAR-index
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 5-00:00:00
#SBATCH -A xxx
#SBATCH -M snowy

ulimit -c unlimited
set -eo pipefail

### Index reference genome in preparation for mapping 
### THIS ONLY NEEDS TO BE DONE ONCE FOR EACH REFERENCE-GENOME/ALIGNER PAIR; 
### REDO IF USING A DIFFERENT VERSION OF STAR ALIGNER
### REDO IF READ LENGTH FOR LIBRARIES OF INTEREST
# load modules
module load bioinfo-tools
module load star/2.7.2b

##### Paths and folders

# Path to folder containing reference genome and annotation file
ref_dir=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/ratGenome_rn6
anno_ref=/crex1/proj/xxx/private/Lili/analysis/bulkRNA/ratGenome_rn6/genes

# reference genome to be indexed
ref_gen=rn6.fa

# Annotation file
anno=rn6.ncbiRefSeq.gtf

#remember initial path
SRCDIR_INI=$(pwd)

############################work with compressed ref and annotation files################
##creat a named pipe(FIFO)
#mkfifo faGz
#zcat ${ref_dir}/${ref_gen} > faGz &

#mkfifo gtfGz
#zcat ${anno_ref}/${anno} > gtfGz &
# index reference genomes

#STAR \
#--runMode genomeGenerate \
#--runThreadN 16 \
#--genomeDir star_index/ \
#--genomeFastaFiles faGz \
#--sjdbGTFfile gtfGz


STAR \
--runMode genomeGenerate \
--runThreadN 10 \
--genomeDir star_index/ \
--genomeFastaFiles ${ref_dir}/${ref_gen} \
--sjdbGTFfile ${anno_ref}/${anno}


#Notes (STAR genomeGenerate)
#--runMode 			genomeGenerate option directs STAR to run genome indices generation job.
#--runThreadN N 	specifies the number of threads that will be used by the program
#--genomeDir 		specifies path to the directory (henceforth called "genome directory" where the
#					genome indices are stored. This directory has to be created (with mkdir) before STAR run
#					and needs to writing permissions. 
#--genomeFastaFiles specified one or more FASTA files with the genome reference sequences
#--sjdbOverhang 	specifies the length of the genomic sequence around the annotated junction
#					to be used in constructing the splice junctions database. Ideally, this length should be equal
#					to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina
#					2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the
#					ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as
#					well as the ideal value. [USE ONLY when annotation file is also provide]
# --sjdbGTFfile 	specifies the path to the file with annotated transcripts in the standard GTF
#					format. STAR will extract splice junctions from this file and use them to greatly improve
#					accuracy of the mapping. While this is optional, and STAR can be run without annotations,
#					using annotations is highly recommended whenever they are available.
# --sjdbOverhang 	specifies the length of the genomic sequence around the annotated junction
#					to be used in constructing the splice junctions database. Ideally, this length should be equal
#					to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina
#					2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the
#					ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as
#					well as the ideal value.

# Annotations in GFF format.
# In addition to the aforementioned options, for GFF3 formatted annotations you need to use
# --sjdbGTFtagExonParentTranscript Parent. In general, for --sjdbGTFfile files STAR only
# processes lines which have --sjdbGTFfeatureExon (=exon by default) in the 3rd field (column). The exons are assigned to the 
# transcripts using parent-child relationship defined by the
# --sjdbGTFtagExonParentTranscript (=transcript id by default) GTF/GFF attribute


# runtime: < 30 min




