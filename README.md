# PolII metrics
This workflow allows calculation of RNA Polymerase II (Pol2) metrics using RNA Pol II ChIP-Seq data:
* Processivity <https://pubmed.ncbi.nlm.nih.gov/32917631/>; 
* Traveling Ratio <https://pubmed.ncbi.nlm.nih.gov/38521065/>;
* Pausing Index <https://pubmed.ncbi.nlm.nih.gov/27259512/>.
The workflow is written in R and includes several Rshiny apps for plotting of calculated metrics and clusters. Please refer to exemplary data files provided in 'data_example_files' folder.

# Rationale
RNA polymerase II (Pol2) transcribes all protein-coding genes, as well as many non-coding genes, and its function is a highly regulated process in mammalian cells. Promoter-proximal pausing of Pol2 after it enters early elongation is an important step in regulating gene expression. The ratio of Pol2 binding over promoter region to the gene body region can estimate the 'pausing' and 'traveling' properties of elongating Pol2. In addition, the Pol2 binding ratio between 5'- and 3'-halves of the gene correlates with processivity of Pol2, e.g. its ability to process the whole gene sequence without 'falling off'. Many elongating factors and therapeutic compounds potentially affect the Pol2 function, and Pol2 ChIP-Seq is a frequently generated datatype to study it. Hence, there is a need to estimate Pol2 metrics under specific experimental conditions. 

# Usage
Typical steps of the workflow include:
* *Optional*: Deriving a subset of actively expressed genes using RNA-Seq counts under your experimental conditions to focus on 'transcriptionally active' Pol 2: *DeriveExpressedGenes()* function. 
	**Note1** If you do not have RNA-Seq data under your experimental conditions, you can use default genes' list ('default_expressed_genes.csv') provided in 'data_example_files' folder. These genes were derived as being expressed in HEK293 cell line. 
	**Note2** *DeriveExpressedGenes()* function expects 4 samples as input: 2 control samples and 2 experimental samples, corresponding to ChIP-Seq design.
* Obtaining genomic coordinates for expressed genes (with length > 2000 bp) from .gff annotation file, promoters and gene body separately, for all three metrics, Processivity (Proc), Traveling Ratio (TR), Pausing Index (PI): *ObtainingBwtoolRegions()* function. This function also saves annotation as variable so we can follow-up on most interesting gene subsets later.
* Calculation of summed ChIP-Seq signal over corresponding genomic coordinates using ChIP-seq data in  .bigWig format and *bwtool* in command line:

bwtool summary PI_denominator_data.bed RPB1_rep1.bigWig output-PI_den_RPB1_rep1.txt  -header -with-sum

    **Note** To install bwtool, follow instructions here <https://github.com/CRG-Barcelona/bwtool>
	
* Importing output .txt files back into R, combining them into single dataframe and normalizing by the corresponding region length using annotation data.
* Using this data as input for *CalculateIndex()* function that outputs: 
	* metric values (average) for every gene in Control vs Experimental condition;
	* significance of statistical comparison Control vs Experimental condition; 
	* quantiles and summary statistics for Control and Experimental condition.
* *Optional*: Lastly, you can cluster genes into 2-4 groups using *kmeansClusteringIndices()* function and follow-up on found gene clusters using *clusters_gene_lengths()* function.

# Built-in Rshiny plots
There are two provided Rshiny plots to visualize your data:
* Overall distribution of a metric of interest (Processivity, Traveling Ratio, Pausing Index) in Control vs Experimental condition;
* Visualization of kmeans-derived gene clusters.