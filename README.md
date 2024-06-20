# PolII_metrics
Calculation of RNA Pol II metrics using Chip-Seq data: Processivity [PMID:32917631], Traveling Ratio [PMID:38521065] and Pausing Index [PMID:27259512].

# Usage

Typical steps include:\
  a) deriving a subset of expressed genes using RNA-Seq counts so to focus on 'transcriptionally engaged' Pol II (optional) (DeriveExpressedGenes() function)\
  b)obtaining genomic regions for these genes (with length > 2000 bp) from .gff annotation file, numerator and denominator separately, for all three indices (ObtainingBwtoolRegions() function). This function also saves annotation file for these genes so we can follow-up on most interesting subsets later\
	c) running bwtool in command line to obtain sums of Chip-Seq signal on corresponding regions:
	(To install bwtool, follow instructions here https://github.com/CRG-Barcelona/bwtool)
	
 bwtool summary PI_denominator_data.bed RPB1_rep1.bigWig output-PI_den_RPB1_rep1.txt  -header -with-sum\
	d) Importing output .txt files into R, combining them into single data.frame and normalizing by region length using annotation data\
	e) Using this data as input for CalculateIndex() function that has 5 outputs: index values (average) for Control vs Treatment; significance of statistical comparison Control vs Treatment; quantiles for both Control and Treatment; visualization\
	f) Lastly, you can try to cluster genes into 2-4 groups (kmeansClusteringIndices() function) and check if there are any differences there.
