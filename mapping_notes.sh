#rat bulk-RNA seq

"""	
	both HISAT and STAR were used to index rn6 genome and align rats bulk-RNA seq
conclusion: 
	STAR performs slightly better than HISAT, STAR will be used in the future
	The article 'Evaluation of Seven Different RNA-Seq Alignment Tools Based on Experimental Data from the Model Plant Arabidopsis thaliana'
	also indicated that HISAT and STAR has similar performance
"""

#########################################################################

"""
	rats genome: mRatBN7 and rn6 were used to STAR analysis
conlusion:
	mRatBN7 and rn6 probably are same since the performance is same,rn6 will be used in the following rats bulk-RNA analysis

"""

#########################################################################

"""
	different annotation files were used in STAR and evaluated by Qualimap
	rn6.ensGene.gtf;
	rn6.ncbiRefSeq.gtf;
	rn6.refGene.gtf; (is more conserved)
	merged_rn6.gtf
conclusion: 
	merged_rn6.gtf is the worst one
	rn6.refGene is the s3econd worst one because of its conserversity
	rn6.ncbiRefSeq is the best one.
	So rn6.ncbiRefSeq.gtf will be used in the future, but please keep alignment based on rn6.ensGene.gtf and rn6.refGene.gtf for the future if wants to check the other traits of genes.
"""

########################################################################

"""
	All above comparisons using both STAR/HISAT and Qualimap summary
"""

########################################################################

"""
where the reads were mapped: rats bulk-RNA seq
	generally low exonic mapping reads, between 55% and 73%, high intronic reads, between 10% and 24%', intergenic 17% -20%
	mapping rates: 74.5% - 84%
	duplicate rates: 11% - 17%

p ossible explanation: 
	1:highintronic mapping may be caused by DNA contamination and/or pre-mRNA
	2:immature transcripts with introns not spliced out yet 2) unannotated exons which might represent alternative splicing
	....
"""






"""
parameters should be considered in assesingmapping quality.
https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html

1: Reads genomic origin: Even if you have high genomic mapping rate for all samples, check to see where the reads are mapping. 
	Ensure that there is not an unusually high number of reads mapping to intronic regions (~30% expected) and fewer than normally observed mapping to exons (~55%). 
	A high intronic mapping suggests possible genomic DNA contamination and/or pre-mRNA.
2: Ribosomal RNA (rRNA) constitutes a large majority of the RNA species in any total RNA preparation. Despite depletion methods, you can never achieve complete rRNA removal. 
	Even with Poly-A enrichment a small percentage of ribosomal RNA can stick to the enrichment beads non-specifically. Excess ribosomal content (> 2%) 
	will normally have to be filtered out so that differences in rRNA mapped reads across samples do not affect alignment rates and skew subsequent normalization of the data.
3: Transcript coverage and 5’-3’ bias: assesing the affect on expression level and on levels of transcript GC content
4: Junction analysis: analysis of junction positions in spliced alignments (i.e known, partly known, novel)
5: Strand specificity: assess the performance of strand-specific library construction methods. 
	The percentage of sense-derived reads is given for each end of the read pair. A non-strand-specific protocol would give values of 50%/50%, 
	whereas strand-specific protocols typically yield 99%/1% or 1%/99% for this metric.
6: the percentage of reads with duplicates (aligning exactly to the same place in the genome) corresponding to  % of reads mapped to multiple loci 
	Usually, duplication levels higher than 20% are not a good sign.
7:  how many (out of all reads) were properly aligned against the reference genome. In the case of bacterial sequencing one would expect >95% successful alignment, 
	but when sequencing a mamallian genome (with many repetitive areas) it may be normal to have as low as 70-80% alignment success. 
	RNA-Seq sequenced regions that are usually well preserved, and thus alignment rates should usually be high.
"""

"""
improve mapping quality:
	1: loose the limit of the minimum reads length, 99 bp was used in the above star mapping
	2: use picard to mark duplicate and check the distribution of those duplicates : up to 50% of duplication can be consider normal to obtain
	3: Number of reads mapped to each chromosome by IdxStats from the Samtools
	4: Gene body coverage, 5'-3'bias
	5: Read distribution across features, exons,introns, intergenic
	6: Estimation of the strandness: Unstranded RNA-Seq data (try three options and then compare the mapped reads)

"""
