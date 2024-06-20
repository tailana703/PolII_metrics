###Calculation of Pol II metrics using Pol II Chip-Seq###
###Processivity_Index [PMID:32917631], Traveling_Ratio [PMID:38521065],
###Pausing_Index [PMID:27259512]

###obtaining regions files from hg38 annotation - will be used for bwtool###

## optional: Filtering expressed genes using RNA-Seq counts and standard pheno-design: 
## 2 Control samples, 2 Treated samples

DeriveExpressedGenes <- function(counts, phenodata) {
  require(dplyr)
  require(DESeq2)
  require(data.table)
  counts_numeric <- mutate_all(counts, function(x) as.numeric(as.character(x)))
  dds <- DESeqDataSetFromMatrix(countData = counts_numeric, colData = pheno, design = ~ group)
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
  dds <- dds[idx,]
  expressed_genes <- rownames(dds)
  expressed_genes <- nth(tstrsplit(expressed_genes, split = "\\."), n =1)
  }

## loading annotation and extracting expressed genes longer than 2000 bp

ObtainingBwtoolRegions <- function(expressed_genes) {
  require(rtracklayer)
  require(dplyr)
  require(data.table)
  ## paste your path to annotation .gtf!!!
  gff_file <- "/.../gencode.v38.primary_assembly.annotation.gtf"
  gff_annot <- rtracklayer::import(gff_file) %>% as.data.frame() 
  gff_annot$gene_id <- nth(tstrsplit(gff_annot$gene_id, split = "\\."), n =1)
                           
  gff_filtered <- gff_annot %>%
    filter(type == "gene") %>%
    filter(gene_id %in% expressed_genes) %>%
    filter(width > 2000)
  
  reduced_gff <- gff_filtered[, c("seqnames", "start", "end", "width", 
                                  "strand", "score", "gene_id")]
  ## Processivity Index
  Proc_numerator_data <- reduced_gff %>%
    mutate(NewStart = start+1000, NewEnd = round(start+width/2))
  Proc_numerator_data <- Proc_numerator_data[,c("seqnames", "NewStart", "NewEnd", "gene_id", "score", "strand")]
  colnames(Proc_numerator_data)[2:3] <- c("start", "end")
  
  Proc_denominator_data <- reduced_gff %>%
    mutate(NewStart = round(start+width/2), NewEnd = end-1000)
  Proc_denominator_data <- Proc_denominator_data[,c("seqnames", "NewStart", "NewEnd", "gene_id", "score", "strand")]
  colnames(Proc_denominator_data)[2:3] <- c("start", "end")
  ## writing as 0-based bed files 
  bedgr1 <- GRanges(Proc_numerator_data[,1], IRanges(Proc_numerator_data[,2]-1, Proc_numerator_data[,3]), Proc_numerator_data$strand)
  export(bedgr1, "Proc_numerator_data.bed", "bed")
  bedgr2 <- GRanges(Proc_denominator_data[,1], IRanges(Proc_denominator_data[,2]-1, Proc_denominator_data[,3]), Proc_denominator_data$strand)
  export(bedgr2, "Proc_denominator_data.bed", "bed")
  
  ## Traveling Ratio
  TR_numerator_data <- reduced_gff %>%
    mutate(NewEnd = start+300)
  TR_numerator_data <- TR_numerator_data[,c("seqnames", "start", "NewEnd", "gene_id", "score", "strand")]
  colnames(TR_numerator_data)[2:3] <- c("start", "end")
  
  TR_denominator_data <- reduced_gff %>%
    mutate(NewStart = start+300)
  TR_denominator_data <- TR_denominator_data[,c("seqnames", "NewStart", "end", "gene_id", "score", "strand")]
  colnames(TR_denominator_data)[2:3] <- c("start", "end")
  
  ## writing as 0-based bed files
  bedgr1 <- GRanges(TR_numerator_data[,1], IRanges(TR_numerator_data[,2]-1, TR_numerator_data[,3]), TR_numerator_data$strand)
  export(bedgr1, "TR_numerator_data.bed", "bed")
  bedgr2 <- GRanges(TR_denominator_data[,1], IRanges(TR_denominator_data[,2]-1, TR_denominator_data[,3]), TR_denominator_data$strand)
  export(bedgr2, "TR_denominator_data.bed", "bed")
  
  ## Pausing Index (PI)
  PI_numerator_data <- reduced_gff %>%
    mutate(NewStart = start-50, NewEnd = start+300)
  PI_numerator_data  <- PI_numerator_data[,c("seqnames", "NewStart", "NewEnd", "gene_id", "score", "strand")]
  colnames(PI_numerator_data )[2:3] <- c("start", "end")
  
  PI_denominator_data <- reduced_gff %>%
    mutate(NewStart = start+300, NewEnd=end +3000)
  PI_denominator_data <- PI_denominator_data[,c("seqnames", "NewStart", "NewEnd", "gene_id", "score", "strand")]
  colnames(PI_denominator_data)[2:3] <- c("start", "end")
  
  ## writing as 0-based bed files
  bedgr1 <- GRanges(PI_numerator_data[,1], IRanges(PI_numerator_data[,2]-1, PI_numerator_data[,3]), PI_numerator_data$strand)
  export(bedgr1, "PI_numerator_data.bed", "bed")
  bedgr2 <- GRanges(PI_denominator_data[,1], IRanges(PI_denominator_data[,2]-1, PI_denominator_data[,3]), PI_denominator_data$strand)
  export(bedgr2, "PI_denominator_data.bed", "bed")

  reduced_gff
}

###Calculating ratios/indices. Default design - unpaired samples (2 controls,
###2 experimental samples), thereby 
### can calculate average values between two biol.replicates

CalculateIndex <- function(df_bwtool_sum, gene, metric) {
  require(reshape2)
  require(ggplot2)
  require(dplyr)
  df_clean <- df_bwtool_sum[complete.cases(df_bwtool_sum), ]
  averages <- data.frame(Group1 = rowMeans(df_clean[,c(1,2)]),
                         Group2 = rowMeans(df_clean[,c(3,4)]))
  averages_clean <- averages[complete.cases(averages),]
  colnames(averages_clean) = c("Control", paste0("sh", gene))
  
  averages_melted <- melt(averages_clean)
  colnames(averages_melted) <- c("Samples", metric)
  
  ## Plotting
  theme_set(theme_bw())
  p <- ggplot(data = averages_melted, aes_string(x="Samples", y=metric, fill = "Samples")) +
    geom_violin() + labs(title="Expressed genes", y=paste0(metric, " ", "Index")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) 
  
  stats <- wilcox.test(as.numeric(averages_clean[,1]), as.numeric(averages_clean[,2]), paired = T, correct = T, na.rm = T)
  res <- list(stats$p.value, quantile(averages_clean[,1], probs = seq(0, 1, 1/4)), 
              quantile(averages_clean[,2], probs = seq(0, 1, 1/4)), p, averages_clean)
  names(res) <- c("p-value", "Stats-Control", paste0("Stats-sh", gene), "plot", "average_indices")
  res
}

###Optional - k-means clustering of indices and comparison of the clusters###

kmeansClusteringIndices <- function(average_indices, ncenters) {
  average_indices[] <- lapply(average_indices, function(x) replace(x, is.infinite(x), NA))
  average_indices <- average_indices[complete.cases(average_indices),]
  cl <- kmeans(average_indices, centers = ncenters, nstart = 10)
  cl
}

################################################################################
###Start
## Importing counts to derive expressed genes
counts <- as.data.frame(read_csv("counts.csv"))
rownames(counts) <- counts$Genes
counts <- counts[, -1]
pheno <- as.data.frame(read_delim("phenodata.txt"))

genes <- DeriveExpressedGenes(counts, pheno)

## Using expressed genes subset, creating regions file corresponding to 
## numerator and denominator for Index calculation (will be in your working directory). 
## Annotation is saved here in order to have it for reference if we want to follow-up on some 
## interesting gene subsets
annot <- ObtainingBwtoolRegions(genes)

## after running bwtool on your bigwig files and obtained regions
## bwtool example (CLI): bwtool summary Proc_numerator_data.bed RPB1_shControl_rep1.bigwig output1.txt  -header -with-sum
## importing results of bwtool here (you have 8 files if there are 2 Control
## samples and 2 Treated samples)

Numerator_shCtrl_rep1 <- read.table("output-PI_num_shControl-rep1.txt", header = FALSE)
Denominator_shCtrl_rep1 <- read.table("output-PI_den_shControl-rep1.txt", header = FALSE)

Numerator_shCtrl_rep2 <- read.table("output-PI_num_shControl-rep2.txt", header = FALSE)
Denominator_shCtrl_rep2 <- read.table("output-PI_den_shControl-rep2.txt", header = FALSE)

Numerator_sh_rep1 <- read.table("output-PI_num_sh-rep1.txt", header = FALSE)
Denominator_sh_rep1 <- read.table("output-PI_den_sh-rep1.txt", header = FALSE)

Numerator_sh_rep2 <- read.table("output-PI_num_sh-rep2.txt", header = FALSE)
Denominator_sh_rep2 <- read.table("output-PI_den_sh-rep2.txt", header = FALSE)

## combining data into single dataframe - Processivity does not require length normalization (5' and 3' regions are of completely the same width),
## but TR and PI do require normalization by numerator or denominator length:
## TR: Length of Numerator is 300, Denominator is annot$width-300 for every gene;
## PI: Length of Numerator is 350, Denominator is annot$width-300+3000 for every gene.
df_bwtool_sum <- data.frame(Control_rep1 = (Numerator_shCtrl_rep1$V10/350)/(Denominator_shCtrl_rep1$V10/(annot$width-300+3000)), 
                            Control_rep2 = (Numerator_shCtrl_rep2$V10/350)/(Denominator_shCtrl_rep2$V10/(annot$width-300+3000)),
                            sh_rep1 = (Numerator_sh_rep1$V10/350)/(Denominator_sh_rep1$V10/(annot$width-300+3000)),
                            sh_rep2 =  (Numerator_sh_rep2$V10/350)/(Denominator_sh_rep2$V10/(annot$width-300+3000)))

## calculate metric you need:
result <-CalculateIndex(df_bwtool_sum, gene = "GOI", 
                        metric = "Pausing_Index")

###explore the results using:
result$average_indices ## Index values
result$`p-value` ## significance of the difference Control vs Treatment 
## Since the power is immense, it is almost always significant.
## check the plot - is the difference noticable at all?
result$`Stats-Control` ## stats on Control Index
result$`Stats-sh` ## stats on Treatment Index
result$plot ## plotting ## feel free to adjust ylim in the ggplot based on 
## your data: scale_y_continuous(limits =  c(XX, XX)) and outlier.shape = NA

###optional - k-means clustering of Indices/genes
## usual number of clusters 2-4
cl <- kmeansClusteringIndices(result$average_indices, 3)
cl$centers ## shows "central points" of the clusters. Are there any differences?
df_ind <-result$average_indices
df_ind[,"cluster"] <- cl$cluster

## Follow-up on the genes within clusters
cl1_PI <- annot[which(df_ind$cluster==1),] ##genes of cluster 1
cl2_PI <- annot[which(df_ind$cluster==2),] ## .. 2
cl3_PI <- annot[which(df_ind$cluster==3),] ## .. 3

## plotting gene lengths
lengths <- data.frame(Cluster = c(rep("Cl1", nrow(cl1_PI)), rep("Cl2", nrow(cl2_PI)), rep("Cl3", nrow(cl3_PI))),
                      Length = c(cl1_PI$width, cl2_PI$width, cl3_PI$width))

p <- ggplot(data =lengths, aes(x=Cluster, y=Length)) +
  geom_boxplot(outlier.shape = NA) + labs(title="Gene lengths breakdown", y="Gene length, bp") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_continuous(limits=c(0, 100000))
p
