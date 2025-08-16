library(GEOquery)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

gse <- getGEO("GSE214647", GSEMatrix = TRUE, AnnotGPL = FALSE)[[1]]
# fetch dataset

sampleInfo <- pData(gse) 

# pData returns a dataframe with sample annotation and metadata

sampleInfo.filtered <- sampleInfo %>% 
  dplyr::select("title", "covid-19 status:ch1") %>%
# extract from sampleInfo the title and covid-19 status only, cuz that's what we'll need for dividing the patients into control and experimental groups
# to select specific columns, we make use of dplyr package, select function.
  dplyr::rename(group = "covid-19 status:ch1",  sampleID = "title")
# we are also renaming the columns that we selected above to group and sampleID

sampleInfo.filtered.ordered <- sampleInfo.filtered %>% dplyr::arrange(sampleInfo.filtered$group)
# we are arranging the rows based on the data in the group (covid-19 status)

sampleInfo.filtered.ordered$SampleType <- factor(rep(c("Covid-negative", "Covid-positive"), c(10, 13)))
# This simply means make the first 10 rows "Covid-negative" and the next 13 "Covid-positive".

sampleInfo.filtered.ordered$SampleType <- relevel(sampleInfo.filtered.ordered$SampleType, "Covid-negative")
# Here we are changing the reference to control. By default, the reference is the one that comes first alphabetically.

sampleInfo.filtered.ordered$SampleAccessions <- rownames(sampleInfo.filtered.ordered)
# Here we are adding sample accessions and setting them = to the row names.

ADData <- read.table("GSE214647_series_matrix.txt", header = T, row.names = 1, check.names = FALSE)
# We are fetching gene expression data
# The file was not loading at first, so I had to delete the weird, unneccessary rows.
# I edited the file such that I only have the genes and IDs.

ADData.ordered <- dplyr::select(ADData, sampleInfo.filtered.ordered$SampleAccessions)
# We are ordering the Gene expression dataframe based on the sampleID from the ordered sampleInfo dataframe.

#########################################################################################
ADdds <- DESeqDataSetFromMatrix(countData = round(ADData.ordered),
                                colData = sampleInfo.filtered.ordered,
                                design = ~ SampleType)
# Creating a DESeq object with the count data and sample info
# I had to convert numbers to integer because I was getting the error: "Some values in assay are not integeres"

View(counts(ADdds)) # Here we are simply viewing the counts to make sure the genes have expression values

ADdds <- DESeq(ADdds)
ADdds.normalized <- estimateSizeFactors(ADdds)
AD.normalized.counts <- counts(ADdds.normalized, normalized=TRUE)
# Running DESeq2 analysis while normalizing the counts.

Contrast <- c("SampleType", "Covid-positive", "Covid-negative")
# Here we are setting the contrast -- that is, Covid-positive vs. Covid-negative

resultsNames(ADdds)
# Shows the columns of the DESeq object (the columns of the differential analysis)

res.table <- results(ADdds, contrast= Contrast, alpha = 0.05)
# Getting the differential analysis results for the specific contrasting groups, based on an alpha of 0.05
summary(res.table)
# Getting a summary of the results --> how many upregulated and how many downregulated, etc.

res.table <- res.table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
# Converting the results to a tidyverse "tibble" with gene names as a column
# Now we have a table with the gene names and their associated statistics
# such as p-val, p-adj, lfc, etc.

padj.cutoff <- 0.05
lfc.cutoff <- 1
# So, we are setting the significance thresholds.

sigOE <- res.table %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
# We are now filtering the differentially expressed genes only

write.csv(sigOE, file="Covid-positive_VS._Covid-negative_DEGs_lfc1_q0.05.csv")
# Then, we are writing the differentially expressed genes to a new Excel file

#########################################################################################
Annotation <- sampleInfo.filtered.ordered %>% 
  dplyr::select(SampleAccessions, SampleType) %>% 
  data.frame(row.names = "SampleAccessions")
# We are making an Annotation dataframe object from a tibble with sampleID and sampleType as columns
# We are making the row names of this dataframe the same as sampleIDs.

vsd <- vst(ADdds, nsub = sum( rowMeans( counts(ADdds, normalized=TRUE)) > 5 ))
# I had to manipulate the nsub because the number of genes is less than 1000, the minimum that is required by vsd.

heat_colors <- brewer.pal(4, "YlOrRd")
# making a vector of 4 colors (based on ylOrRd palette) to be used in the heatmap
DEGs <- sigOE[, 1]
#Extracting all the gene names from the differentially expressed genes as a dataframe

DEGs <- DEGs[["gene"]]
# Extracting the individal gene names themselves

PH = pheatmap(assay(vsd)[DEGs, ],
         color = heat_colors, 
         cluster_rows = T, 
         cluster_cols = F,
         show_rownames = T,
         annotation = Annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 6, 
         height = 400)

# Now we are drawing a heat map with the specified parameters

EV = EnhancedVolcano(res.table, 
                lab = "", 
                x = 'log2FoldChange',
                y = 'padj',
                title = "Covid-positive Vs. Covid-negative",
                pCutoff = 0.05,
                FCcutoff = 1.0)
# Now we are drawing a volcano plot with the specified paramters

ggsave("EnhancedVolcano.jpg", plot=EV)
ggsave("HeatMap.jpg", plot=PH)
