###Calculation of Pol II metrics using Pol II Chip-Seq###
###Processivity_Index [PMID:32917631], Traveling_Ratio [PMID:38521065],
###Pausing_Index [PMID:27259512]

###Step1: obtaining region files from gtf annotation - will be used for bwtool###

##Step1.1 (optional): Filtering expressed genes using RNA-Seq counts and standard design: 
## 2 Control samples, 2 Experimental samples

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

##Step1.2: loading annotation and extracting unique expressed genes longer than 2000 bp
##Make sure to paste path to annotation file
ObtainingBwtoolRegions <- function(expressed_genes, gff_file) {
  require(rtracklayer)
  require(dplyr)
  require(data.table)
  gff_annot <- rtracklayer::import(gff_file) %>% as.data.frame() 
  gff_annot$gene_id <- nth(tstrsplit(gff_annot$gene_id, split = "\\."), n =1)
                           
  gff_filtered <- gff_annot %>%
    filter(type == "gene") %>%
    filter(gene_id %in% expressed_genes) %>%
    filter(!duplicated(gene_id)) %>%
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

###Step 2: Calculating ratios/indices. Default design - unpaired samples (2 controls,
###2 experimental samples), thereby we can calculate average values between two biol.replicates

CalculateIndex <- function(df_bwtool_sum, gene) {
  require(reshape2)
  require(ggplot2)
  require(dplyr)
  df_clean <- df_bwtool_sum[complete.cases(df_bwtool_sum), ]
  averages <- data.frame(Group1 = rowMeans(df_clean[,c(1,2)]),
                         Group2 = rowMeans(df_clean[,c(3,4)]))
  averages_clean <- averages[complete.cases(averages),]
  averages_clean <- averages_clean[!is.infinite(averages_clean[,1]),]
  averages_clean <- averages_clean[!is.infinite(averages_clean[,2]),]
  colnames(averages_clean) = c("Control", paste0("sh", gene))
  
  stats <- wilcox.test(as.numeric(averages_clean[,1]), as.numeric(averages_clean[,2]), paired = T, correct = T, na.rm = T)
  res <- list(stats$p.value, quantile(averages_clean[,1], probs = seq(0, 1, 1/4)), 
              quantile(averages_clean[,2], probs = seq(0, 1, 1/4)), averages_clean)
  names(res) <- c("p-value", "Stats-Control", paste0("Stats-sh", gene), "average_indices")
  res
}

###Step 3 (optional) - k-means clustering of indices and comparison of the clusters###

kmeansClusteringIndices <- function(average_indices, ncenters) {
  average_indices[] <- lapply(average_indices, function(x) replace(x, is.infinite(x), NA))
  average_indices <- average_indices[complete.cases(average_indices),]
  cl <- kmeans(average_indices, centers = ncenters, nstart = 10)
  cl
}

###Step4 (optional) - follow-up on found gene clusters: gene lengths 

clusters_gene_lengths <- function(clusters = cl, annotation = annot) {
  cluster_genes_names <- list()
  cluster_genes_lengths <- list()
  for (i in 1:length(cl$size)) {
    cluster_genes_names[[i]] <- names(cl$cluster[cl$cluster==i])
    nrow_annot <- which(annot$gene_id %in% names(cl$cluster[cl$cluster==i]))
    cluster_genes_lengths[[i]] <- annot$width[nrow_annot]
  }
  result <- list(cluster_genes_names, cluster_genes_lengths)
  names(result) <- c("genes", "lengths")
  result
}


################################################################################
###Start

#setwd("//") ##specify your working directory
# In case of difficulties when installing packages, use the following code (example - installing DESEq2 package):
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("DESeq2")

## Importing counts to derive expressed genes, example files can be found in the repo
library(readr)
counts <- as.data.frame(read_csv("counts.csv"))

## if counts are stored in .xlsx or .txt format, use different functions:
#library(readxl)
#library(readr)
#counts <- as.data.frame(read_excel("counts.xlsx"))
#counts <- as.data.frame(read_delim("counts.txt"))

rownames(counts) <- counts$Genes
counts <- counts[, -1]
pheno <- as.data.frame(read_delim("phenodata.txt"))
pheno$group <- as.factor(pheno$group)
genes <- DeriveExpressedGenes(counts, pheno)
head(genes)

##alternatively, use default gene set that is provided in the repo
default_genes <- read_csv("default_expressed_genes.csv") 
genes <- default_genes$genes

## Using expressed genes subset, creating .bed regions file corresponding to 
## numerator and denominator for Index calculation (will be in your working directory). 
## Make sure .gtf annotation is on your hard drive and specify the path:

#gff_file <- "/Users/sveta/gencode.v19.annotation.gtf" #OR
gff_file <- "/Users/sveta/gencode.v38.primary_assembly.annotation.gtf"
annot <- ObtainingBwtoolRegions(genes, gff_file)

## after running bwtool on your bigwig files and obtained regions
## bwtool example (CLI): bwtool summary Proc_numerator_data.bed RPB1_shControl_rep1.bigwig output1.txt  -header -with-sum
## importing results of bwtool here (you have 8 files per Index if there are 2 Control
## samples and 2 Treated samples)

PI_Numerator_Ctrl_rep1 <- read.table("output-PI_num_shControl-rep1.txt",header = FALSE)
PI_Denominator_Ctrl_rep1 <- read.table("output-PI_den_shControl-rep1.txt", header = FALSE)
PI_Numerator_Ctrl_rep2 <- read.table("output-PI_num_shControl-rep2.txt", header = FALSE)
PI_Denominator_Ctrl_rep2 <- read.table("output-PI_den_shControl-rep2.txt", header = FALSE)
PI_Numerator_Exp_rep1 <- read.table("output-PI_num_shGOI-rep1.txt", header = FALSE)
PI_Denominator_Exp_rep1 <- read.table("output-PI_den_shGOI-rep1.txt", header = FALSE)
PI_Numerator_Exp_rep2 <- read.table("output-PI_num_shGOI-rep2.txt", header = FALSE)
PI_Denominator_Exp_rep2 <- read.table("output-PI_den_shGOI-rep2.txt", header = FALSE)

TR_Numerator_Ctrl_rep1 <- read.table("output-TR_num_shControl-rep1.txt", header = FALSE)
TR_Denominator_Ctrl_rep1 <- read.table("output-TR_den_shControl-rep1.txt", header = FALSE)
TR_Numerator_Ctrl_rep2 <- read.table("output-TR_num_shControl-rep2.txt", header = FALSE)
TR_Denominator_Ctrl_rep2 <- read.table("output-TR_den_shControl-rep2.txt", header = FALSE)
TR_Numerator_Exp_rep1 <- read.table("output-TR_num_shGOI-rep1.txt", header = FALSE)
TR_Denominator_Exp_rep1 <- read.table("output-TR_den_shGOI-rep1.txt", header = FALSE)
TR_Numerator_Exp_rep2 <- read.table("output-TR_num_shGOI-rep2.txt", header = FALSE)
TR_Denominator_Exp_rep2 <- read.table("output-TR_den_shGOI-rep2.txt", header = FALSE)

Proc_Numerator_Ctrl_rep1 <- read.table("output-Proc_num_shControl-rep1.txt", header = FALSE)
Proc_Denominator_Ctrl_rep1 <- read.table("output-Proc_den_shControl-rep1.txt", header = FALSE)
Proc_Numerator_Ctrl_rep2 <- read.table("output-Proc_num_shControl-rep2.txt", header = FALSE)
Proc_Denominator_Ctrl_rep2 <- read.table("output-Proc_den_shControl-rep2.txt", header = FALSE)
Proc_Numerator_Exp_rep1 <- read.table("output-Proc_num_shGOI-rep1.txt", header = FALSE)
Proc_Denominator_Exp_rep1 <- read.table("output-Proc_den_shGOI-rep1.txt", header = FALSE)
Proc_Numerator_Exp_rep2 <- read.table("output-Proc_num_shGOI-rep2.txt", header = FALSE)
Proc_Denominator_Exp_rep2 <- read.table("output-Proc_den_shGOI-rep2.txt", header = FALSE)

## combining data into three dataframes corresponding to Index and normalization by region lengths
df_bwtool_PI <- data.frame(Control_rep1 = (PI_Numerator_Ctrl_rep1$V10/350)/(PI_Denominator_Ctrl_rep1$V10/(annot$width-300+3000)), 
                            Control_rep2 = (PI_Numerator_Ctrl_rep2$V10/350)/(PI_Denominator_Ctrl_rep2$V10/(annot$width-300+3000)),
                            Exp_rep1 = (PI_Numerator_Exp_rep1$V10/350)/(PI_Denominator_Exp_rep1$V10/(annot$width-300+3000)),
                            Exp_rep2 =  (PI_Numerator_Exp_rep2$V10/350)/(PI_Denominator_Exp_rep2$V10/(annot$width-300+3000)))
df_bwtool_TR <- data.frame(Control_rep1 = (TR_Numerator_Ctrl_rep1$V10/300)/(TR_Denominator_Ctrl_rep1$V10/(annot$width-300)), 
                           Control_rep2 = (TR_Numerator_Ctrl_rep2$V10/300)/(TR_Denominator_Ctrl_rep2$V10/(annot$width-300)),
                           Exp_rep1 = (TR_Numerator_Exp_rep1$V10/300)/(TR_Denominator_Exp_rep1$V10/(annot$width-300)),
                           Exp_rep2 =  (TR_Numerator_Exp_rep2$V10/300)/(TR_Denominator_Exp_rep2$V10/(annot$width-300)))
df_bwtool_Proc <- data.frame(Control_rep1 = (Proc_Numerator_Ctrl_rep1$V10)/(Proc_Denominator_Ctrl_rep1$V10), 
                           Control_rep2 = (Proc_Numerator_Ctrl_rep2$V10)/(Proc_Denominator_Ctrl_rep2$V10),
                           Exp_rep1 = (Proc_Numerator_Exp_rep1$V10)/(Proc_Denominator_Exp_rep1$V10),
                           Exp_rep2 =  (Proc_Numerator_Exp_rep2$V10)/(Proc_Denominator_Exp_rep2$V10))

rownames(df_bwtool_PI) <- annot$gene_id
rownames(df_bwtool_TR) <- annot$gene_id
rownames(df_bwtool_Proc) <- annot$gene_id

## calculate metric (do not forget to change the name of the gene):
gene <- "GOI" ## paste the name of your knockdown gene (GOI = gene of interest) or treatment here
PI <-CalculateIndex(df_bwtool_PI, gene = gene)
TR <-CalculateIndex(df_bwtool_TR, gene = gene)
Proc <-CalculateIndex(df_bwtool_Proc, gene = gene)

###explore the results using:
head(PI$average_indices) 
PI$`p-value` ## significance of the difference Control vs Treatment 
## Since the power is immense, it is almost always significant.
PI$`Stats-Control` ## stats on Control Index
PI$`Stats-shGOI` ## stats on Treatment/Experimental Index

##Rshiny app to plot average indices (can customize plotting parameters and which index to show)##
library(shiny)
# Define UI for application that draws a boxplot
ui <- fluidPage(
  titlePanel("PolII metrics"),
  sidebarLayout(
    sidebarPanel(selectInput("metric", "Metric:", choices = c("Pausing Index", 
                                                              "Traveling Ratio",
                                                              "Processivity Index")),
                 selectInput("color", "Color:", choices = c("white", "orange", "darkgreen", "blue", "red")),
                 checkboxInput("outliers", "Show outliers", TRUE)
    ),
    mainPanel(
      plotOutput("boxPlot")
    )
  )
)


# Define server logic required to draw a boxplot
server <- function(input, output, session) {
  x <- reactive({input$metric})

  output$boxPlot <- renderPlot({
    if (x() == "Pausing Index") {
      values <- PI$average_indices
    }
    if (x() == "Traveling Ratio") {
      values <- TR$average_indices
    }
    if (x() == "Processivity Index") {
      values <- Proc$average_indices
    }
    boxplot(values,
            outline = input$outliers,
            col = input$color, pch = 19,
            main = paste("PolII metric:", input$metric),
            xlab = input$metric,
            horizontal = T,
            boxwex = 0.25)
  })
  
}

shinyApp(ui = ui, server = server)

###optional - k-means clustering of Indices/genes
## usual number of clusters 2-4 - put the desired number in the call

cl <- kmeansClusteringIndices(PI$average_indices, ncenters = 3)
cl$centers ## shows "central points" of the clusters. Are there any differences?

## Follow-up on the genes within clusters using clusters_gene_lengths function
resulting_lengths <- clusters_gene_lengths(clusters = cl, annotation = annot)

## example of plotting gene lengths when number of clusters = 3:
boxplot(resulting_lengths$lengths[[1]], 
        resulting_lengths$lengths[[2]], 
        resulting_lengths$lengths[[3]], 
        outline = F, col = "purple",
        ylab = "Gene length, bp",
        xlab = "Clusters",
        main = "PI cluster vs gene length",
        names = c("1", "2", "3"))


###Bonus: cluster visualization with Rshiny!
library(rmdexamples)
kmeans_cluster(PI$average_indices)
kmeans_cluster(TR$average_indices)
kmeans_cluster(Proc$average_indices)

