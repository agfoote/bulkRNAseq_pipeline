#required packages for analysis 
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)
library(pheatmap)
library(EnhancedVolcano)
library(topGO)
library(apeglm)
setwd("/Users/alexanderfoote/R_Projects/bulkRNA_analysis_IAV_21dpi")

# load data and remove extra header stuff. also version numbers 
counts <-read.table("/Users/alexanderfoote/R_Projects/bulkRNA_analysis_IAV_21dpi/featurecounts_IAV_21dpi.txt",header =T)
counts$Geneid <- sub("\\.\\d+", "",counts$Geneid)
#change row names to gene ID
row.names(counts) <- counts$Geneid
#exclude extra columns
counts <- counts[,-c(1:6)]
#rename column names 
colnames(counts) <- c("Saline_1","Saline_2","Saline_3","Saline_4","IAV_1","IAV_2","IAV_3","IAV_4","IAV_5","IAV_6","IAV_7")
# Remove the IAV_4 column since this seems to be an outlier
counts <- counts[, colnames(counts) != "IAV_4"]
colnames(counts)
#load conditions 
coldata<- read.csv("/Users/alexanderfoote/R_Projects/bulkRNA_analysis_IAV_21dpi/RNAseq_coldata_IAV_21dpi_v2.csv",row.names = 1, header=TRUE)
rownames(coldata)
#converting to factors for Deseq2 input 
coldata$condition<-factor(coldata$condition)
#check to see if rownames of coldata matches colnames of counts. should yield a line of TRUES 
colnames(counts) == rownames(coldata)


#Deseq2 Analysis
dds <-DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <-results(dds)
summary(res)

#Adding gene names 
res$symbol <- mapIds(org.Mm.eg.db, keys= row.names(counts), column= "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res_ordered <- res[order(res$pvalue),]

#retrieving significant genes 
head(res_ordered)
resSig05 <- subset(as.data.frame(res), padj < 0.05)
resSig01 <- subset(as.data.frame(res), padj <0.1)
#reording by padj
resSig05 <-resSig05[order(resSig05$padj),]
resSig01 <-resSig01[order(resSig01$padj),]
#in total 
setwd("/Users/alexanderfoote/R_Projects/bulkRNA_analysis_IAV_21dpi/csv")
nrow(resSig05)
res_05_list<-write.csv(resSig05,"IAV_21dpi_v2.csv")

#cleaning up the data
resSig05<-resSig05[complete.cases(resSig05),]

#optional- viewing counts of a specific ensembl gene 
plotCounts(dds, gene="ENSMUSG00000021732", intgroup="condition")

#PCA- global view of the data 
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)

data <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Plot PCA")

#labeling samples 
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Plot PCA") +
  geom_label_repel(aes(label=name))

#plotting MA counts, log2 fold changes over mean of normalized counts 
DESeq2::plotMA(res, ylim=c(-2,2), main="DESeq2: Control&IAV plot")
summary(res)

## check to see number of low counts that might contribute to strangeness 
## applying a lfc shrink 
resultsNames(dds)[2]
resLFC <- lfcShrink(dds, coef="condition_IAV_vs_Control", type="apeglm")
DESeq2::plotMA(resLFC, ylim=c(-2,2), main="IAV plot: LFC Shrink" )

resLFC$symbol <- mapIds(org.Mm.eg.db, keys= row.names(counts), column= "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res_LFC_ordered <- resLFC[order(resLFC$pvalue),]

#pheatmap 
#take top 20 
top_genes<- resSig05[1:20,]
mat<- assay(rld)[row.names(top_genes),]
row.names(mat) <- top_genes$symbol
df <- as.data.frame(coldata)
setwd("/Users/alexanderfoote/R_Projects/bulkRNA_analysis_IAV_21dpi/plots")
name <- "IAV_21dpi_v2"
pdf(paste(name,"_heatmap.pdf",sep=""))
pheatmap(mat, annotation_col=df,scale = "row")
dev.off()

#sample to sample distance 
sampleDist <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDist)

rownames(sampleDistMatrix)<- paste(rld$condition)
colnames(sampleDistMatrix)<- paste(rld$condition)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, cluster_distance_rows= sampleDist, clustering_distance_cols = sampleDist, col=colors)

# volcano plot  
setwd("/Users/alexanderfoote/R_Projects/bulkRNA_analysis_IAV_21dpi/plots")
name <- "IAV_21dpi_v2"
pdf(paste(name,"_Volcanoplot.pdf",sep=""))
EnhancedVolcano(res, lab=res$symbol,x='log2FoldChange', 
                y='pvalue', 
                FCcutoff = 0.2, 
                pCutoff = 10e-4,
                selectLab = c("H2-Ea",'Gzmb','Gzma','Ccl5','Cd8a','Cd8b1','Myl4','Lep','Sult1b1','Prf1',"Mcpt1","Sftpa1","Sult1c1","Cckar","Tacr1"),
                title='IAV Cell Markers - 21dpi',
                subtitle= 'Differential Expression', 
                caption= bquote(~Log[2]~ "fold change cutoff, 0.2; p-value cutoff, 0.05"),
                xlim = c(-5,5),
                ylim = c(0,35),
                pointSize = 1.5,
                labSize = 5,
                max.overlaps = 35,
                drawConnectors = TRUE,
                endsConnectors = 'last',
                col=c('grey','grey','grey','red'),
                legendPosition = 'right',
                legendLabSize = 12,
                widthConnectors = 0.5)
dev.off()