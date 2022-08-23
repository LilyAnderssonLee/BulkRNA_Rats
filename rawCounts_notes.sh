"""
featurecounts: generate tables of gene counts; It ignores reads mapping equally well to multiple positions
	1: unstrand-specific reads;
	2: strand-specific reads;
	3: reverse strand-specific reads

conclusion: unstrand-specific performs best based on summary comparisons. So we used unstranded library. 
"""

#############################################################################

"""
htseq:
	conduct on each bam file and then combine them
"""

