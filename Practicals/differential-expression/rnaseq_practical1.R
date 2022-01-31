## Install packages from CRAN
install.packages("ggsci")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("ggrepel")
install.packages("RColorBrewer")
install.packages("pheatmap")
install.packages("VennDiagram")
install.packages("gridExtra")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("biomaRt")

# Define the base directory.
base.dir <- ""

# Create the new directory
dir.create(paste(base.dir, "RNASeqPracticalSesion", sep=""))

# Change the working directory
setwd(paste(base.dir, "RNASeqPracticalSesion", sep=""))
dir.create("Data")

# Get all file names with a certain pattern from a directory
fileList <- list.files("./Data/GSE69549_RAW/", pattern = "", full.names = TRUE)
# For each of the components of the fileList, read them into R, and arrange them in a list.
countList <- lapply(fileList, read.delim, sep="\t", header=TRUE, row.names=1)
# Take each the loaded matrix and merge the by it's column (cbind)
countTable <- do.call(cbind, countList)

# Check for the dimensions of the count table. ROWS x COLUMS
dim(countTable)

# Print the firs 5 rows, by the first 5 columns
countTable[1:5, 1:5]

# remove the countList object
rm(countList)

# Print the"head" of the column names.
head(colnames(countTable))

# Use the information in the column names to generate a data.frame comprising the sample information.
colData <- data.frame(row.names = colnames(countTable),
                      type= unlist(lapply(colnames(countTable),
                                          function(x){unlist(strsplit(x, split = "_"))[1]})),
                      sample= unlist(lapply(colnames(countTable),
                                            function(x){unlist(strsplit(x, split = "_"))[2]})),
                      stimulation= unlist(lapply(colnames(countTable),
                                                 function(x){unlist(strsplit(x, split = "_"))[3]}))
)

# Add library size, which is the total amount of gene reads per sample
colData$libSize <- colSums(countTable)
# First 10 rows of the colData object
head(colData)

# Make sure that the colnames of your count table are matching with the row names of your colData
all(rownames(colData) %in% colnames(countTable))

# Save your countTable and colData
write.csv(countTable, file="./Data/countTable.csv")
write.csv(colData, file="./Data/colData.csv")

# table generates a table of frequencies given a vector, we can check how many samples we should have per stimulation
table(colData$stimulation)

# It's possible to include two vectors and count frequencies given two conditions.
table(colData$stimulation,colData$type)

unstim_sampleNames <- rownames(colData)[which(colData$stimulation == "UNS")]
unstim_countTable <- countTable[,unstim_sampleNames]
unstim_colData <- colData[unstim_sampleNames,] #colData has sample names by row

cd3_sampleNames <- rownames(colData)[which(colData$stimulation == "CD3")]
cd3_countTable <- countTable[,cd3_sampleNames]
cd3_colData <- colData[cd3_sampleNames,] #colData has sample names by row

pma_sampleNames <- rownames(colData)[which(colData$stimulation == "PMA")]
pma_countTable <- countTable[,pma_sampleNames]
pma_colData <- colData[pma_sampleNames,] #colData has sample names by row

# Load packages
library(DESeq2)
library(ggplot2)
library(ggsci)
library(ggrepel)

dds_unstim <- DESeqDataSetFromMatrix(countData = unstim_countTable,
                              colData = unstim_colData,
                              design = ~ type+sample)

dds_cd3 <- DESeqDataSetFromMatrix(countData = cd3_countTable,
                              colData = cd3_colData,
                              design = ~ type+sample)

dds_pma <- DESeqDataSetFromMatrix(countData = pma_countTable,
                              colData = pma_colData,
                              design = ~ type+sample)


vst_unstim <- assay(vst(dds_unstim, blind=FALSE))
vst_cd3 <- assay(vst(dds_cd3, blind=FALSE))
vst_pma <- assay(vst(dds_pma, blind=FALSE))

# To calculate the components by sample we need to transpose our matrix of normalized gene expression
pcData <- prcomp(t(vst_unstim))
pcVar <- summary(pcData)
# By getting the summary() of a prcomp object (in this case pcData) we can also obtain the total amount of variance explained by each of the components.
pcVar$importance

# We can then extract the variance explained by components 1 and 2.
varPC1 <- pcVar$importance[2,1]
varPC2 <- pcVar$importance[2,2]

pcPlotData <- data.frame(pcData$x[,1:4], colData[rownames(pcData$x),])
pcaPlot_unstim <- ggplot(pcPlotData, aes(x=PC1 , y=PC2 , color=type))+
                  geom_jitter(alpha=0.6)+
                  facet_grid(~stimulation, scales = "free")+
                  xlab(paste0("PC1 explained variance = ", varPC1*100, "%"))+
                  ylab(paste0("PC2 explained variance = ", varPC2*100, "%"))+
                  scale_color_aaas()+
                  theme_bw()+
                  theme(legend.position = "bottom")+
                  guides(col = guide_legend(ncol = 8))


pcData <- prcomp(t(vst_cd3))
pcVar <- summary(pcData)
varPC1 <- pcVar$importance[2,1]
varPC2 <- pcVar$importance[2,2]

pcPlotData <- data.frame(pcData$x[,1:4], colData[rownames(pcData$x),])
pcaPlot_cd3 <- ggplot(pcPlotData, aes(x=PC1 , y=PC2 , color=type))+
                  geom_jitter(alpha=0.6)+
                  facet_grid(~stimulation, scales = "free")+
                  xlab(paste0("PC1 explained variance = ", varPC1*100, "%"))+
                  ylab(paste0("PC2 explained variance = ", varPC2*100, "%"))+
                  scale_color_aaas()+
                  theme_bw()+
                  theme(legend.position = "bottom")+
                  guides(col = guide_legend(ncol = 8))


pcData <- prcomp(t(vst_pma))
pcVar <- summary(pcData)
varPC1 <- pcVar$importance[2,1]
varPC2 <- pcVar$importance[2,2]

pcPlotData <- data.frame(pcData$x[,1:4], colData[rownames(pcData$x),])
pcaPlot_pma <- ggplot(pcPlotData, aes(x=PC1 , y=PC2 , color=type))+
                  geom_jitter(alpha=0.6)+
                  facet_grid(~stimulation, scales = "free")+
                  xlab(paste0("PC1 explained variance = ", varPC1*100, "%"))+
                  ylab(paste0("PC2 explained variance = ", varPC2*100, "%"))+
                  scale_color_aaas()+
                  theme_bw()+
                  theme(legend.position = "bottom")+
                  guides(col = guide_legend(ncol = 8))

library(gridExtra)
grid.arrange(pcaPlot_unstim, pcaPlot_cd3, pcaPlot_pma, nrow=1)

library(RColorBrewer)
library(pheatmap)

# Again we need to transpose our matrix to then calculate the distance between each of the samples.
sampleDists <- dist(t(vst_unstim))
sampleDistMatrix <- as.matrix(sampleDists)

# By using brewer.pal() we can generate a palette of colors, for more colors check (http://colorbrewer2.org/)
colors <- colorRampPalette(brewer.pal(9, "GnBu"))(255)

pheatmap(sampleDistMatrix, main = "Unstimulated samples",
         show_colnames = FALSE,
         annotation = unstim_colData[,c("stimulation","type")],
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

########## Heatmap for CD3 samples
sampleDists <- dist(t(vst_cd3))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix, main = "CD3 stimulated samples",
         show_colnames = FALSE,
         annotation = cd3_colData[,c("stimulation","type")],
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

########## Heatmap for PNA samples
sampleDists <- dist(t(vst_pma))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix, main = "PMA stimulated samples",
         show_colnames = FALSE,
         annotation = pma_colData[,c("stimulation","type")],
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Load
library(DESeq2)
library(ggplot2)
library(ggsci)
library(ggrepel)

# To properly compare control versus celiac, we need to define that in our colData objects, therefore we need to set the colData$type as a factor, and the firs level should be control

unstim_colData$type <- factor(unstim_colData$type, levels= c("Control", "Coeliac"))
cd3_colData$type <- factor(cd3_colData$type, levels= c("Control", "Coeliac"))
pma_colData$type <- factor(pma_colData$type, levels= c("Control", "Coeliac"))

# We now generate the new DESeq objects and run the differential expression analysis.
dds_unstim <- DESeqDataSetFromMatrix(countData = unstim_countTable,
                              colData = unstim_colData,
                              design = ~ type)

dds_unstim <- DESeq(dds_unstim)
res_unstim <- results(dds_unstim)
head(res_unstim)

# Given the conditions we declared for differential expression analysis we can subset our list of
de_unstim <- res_unstim[which(res_unstim$padj <= 0.05),]
nrow(de_unstim)

dds_cd3 <- DESeqDataSetFromMatrix(countData = cd3_countTable,
                              colData = cd3_colData,
                              design = ~ type)
dds_cd3 <- DESeq(dds_cd3)
res_cd3 <- results(dds_cd3)
de_cd3 <- res_unstim [which(res_cd3$padj <= 0.05),]
nrow(de_cd3)


dds_pma <- DESeqDataSetFromMatrix(countData = pma_countTable,
                              colData = pma_colData,
                              design = ~ type)
dds_pma <- DESeq(dds_pma)
res_pma <- results(dds_pma)
de_pma <- res_unstim [which(res_pma$padj <= 0.05),]
nrow(de_pma)

library(VennDiagram)

DEG_list <- list(Unstimulated= rownames(de_unstim),
                 CD3=rownames(de_cd3),
                 PMA=rownames(de_pma)
                 )
vennDiagram <- venn.diagram(DEG_list,
                            filename = NULL,
                            main = "Shared Differential expressed genes between treatments"
                            )
plot.new()
grid.draw(vennDiagram)

# Lets create a new directory to save our results
dir.create("./Results/")

# The write.csv function will generate an Excel "friendly" file.
unstim_fileName <- "./Results/diffExpGenes_unstim.csv"
write.csv(de_unstim, unstim_fileName)

cd3_fileName <- "./Results/diffExpGenes_cd3.csv"
write.csv(de_cd3, cd3_fileName)

pma_fileName <- "./Results/diffExpGenes_pma.csv"
write.csv(de_pma, pma_fileName)

library(gridExtra)

pData <- as.data.frame(res_unstim[which(!is.na(res_unstim$padj)),])
unstim_Volcano <- ggplot(pData, aes(x=log2FoldChange, y= -log10(padj)))+
                    geom_point(aes(color= padj <= 0.05))+
                    geom_hline(yintercept = 0, lwd=1, alpha= 0.6)+
                    geom_vline(xintercept = 0, lwd=1, alpha= 0.6)+
                    scale_color_d3()+
                    ggtitle("Unstimulated")+
                    theme_bw()

print(unstim_Volcano)

pData$top10label <- NA
pData$top10label[order(pData$padj)[1:10]] <- rownames(pData)[order(pData$padj)[1:10]]
unstim_Volcano <- ggplot(pData, aes(x=log2FoldChange, y= -log10(padj)))+
                    geom_point(aes(color= padj <= 0.05))+
                    geom_hline(yintercept = 0, lwd=1, alpha= 0.6)+
                    geom_vline(xintercept = 0, lwd=1, alpha= 0.6)+
                    scale_color_d3()+
                    ggtitle("Unstimulated")+
                    geom_text_repel(aes(label=top10label))+ ##add the lables in the top 10
                    theme_bw()+
                    theme(legend.position = "bottom")

pData <- as.data.frame(res_cd3[which(!is.na(res_cd3$padj)),])
pData$top10label <- NA
pData$top10label[order(pData$padj)[1:10]] <- rownames(pData)[order(pData$padj)[1:10]]
cd3_Volcano <- ggplot(pData, aes(x=log2FoldChange, y= -log10(padj)))+
                    geom_point(aes(color= padj <= 0.05))+
                    geom_hline(yintercept = 0, lwd=1, alpha= 0.6)+
                    geom_vline(xintercept = 0, lwd=1, alpha= 0.6)+
                    scale_color_jama()+
                    geom_text_repel(aes(label=top10label))+ ##add the lables in the top 10
                    ggtitle("CD3 stimulated")+
                    theme_bw()+
                    theme(legend.position = "bottom")


pData <- as.data.frame(res_pma[which(!is.na(res_pma$padj)),])
pData$top10label <- NA
pData$top10label[order(pData$padj)[1:10]] <- rownames(pData)[order(pData$padj)[1:10]]
pma_Volcano <- ggplot(pData, aes(x=log2FoldChange, y= -log10(padj)))+
                    geom_point(aes(color= padj <= 0.05))+
                    geom_hline(yintercept = 0, lwd=1, alpha= 0.6)+
                    geom_vline(xintercept = 0, lwd=1, alpha= 0.6)+
                    scale_color_rickandmorty()+
                    geom_text_repel(aes(label=top10label))+ ##add the lables in the top 10
                    ggtitle("PMA stimulated")+
                    theme_bw()+
                    theme(legend.position = "bottom")

grid.arrange(unstim_Volcano, cd3_Volcano, pma_Volcano, nrow=1)

library(pheatmap)

vst_unstim <- assay(vst(dds_unstim))
deGenes <- vst_unstim[rownames(de_unstim),]
pheatmap(deGenes, scale = "row",show_rownames = FALSE, main = "Differential expressed genes unstimulated")

vst_cd3 <- assay(vst(dds_cd3))
deGenes <- vst_cd3[rownames(de_cd3),]
pheatmap(deGenes, scale = "row",show_rownames = FALSE, main = "Differential expressed genes CD3")

vst_pma <- assay(vst(dds_pma))
deGenes <- vst_pma[rownames(de_pma),]
pheatmap(deGenes, scale = "row",show_rownames = FALSE, main = "Differential expressed genes PMA")

library(reshape2)

top10Genes <- rownames(de_unstim)[order(de_unstim$padj)[1:10]]
# Using melt to literaly melt a wide data.frame into a long data.frame
pData <- melt(vst_unstim[top10Genes,])
# Add your sample information.
pData <- cbind(pData, unstim_colData[as.character(pData$Var2),])

top10_Plot <- ggplot(pData, aes(x= type, y= value))+
                geom_jitter(alpha=0.8)+
                geom_boxplot(alpha=0.6)+
                facet_grid(~Var1, scale="free")+
                ylab("VST expression values")+
                xlab("")+
                theme_bw()+
                theme(axis.text.x = element_text(angle=45, hjust = 1))
print(top10_Plot)

top10Genes <- rownames(de_cd3)[order(de_cd3$padj)[1:10]]
# Using melt to literaly melt a wide data.frame into a long data.frame
pData <- melt(vst_cd3[top10Genes,])
# Add your sample information.
pData <- cbind(pData, cd3_colData[as.character(pData$Var2),])

top10_Plot <- ggplot(pData, aes(x= type, y= value))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(alpha=0.6)+
  facet_grid(~Var1, scale="free")+
  ylab("CD3 expression values")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))
print(top10_Plot)

top10Genes <- rownames(de_pma)[order(de_pma$padj)[1:10]]
# Using melt to literaly melt a wide data.frame into a long data.frame
pData <- melt(vst_pma[top10Genes,])
# Add your sample information.
pData <- cbind(pData, pma_colData[as.character(pData$Var2),])

top10_Plot <- ggplot(pData, aes(x= type, y= value))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(alpha=0.6)+
  facet_grid(~Var1, scale="free")+
  ylab("PMA expression values")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))
print(top10_Plot)

# The write.csv function will generate an Excel "friendly" file.
unstim_fileName <- "./Results/diffExpGenes_unstim_all.csv"
write.csv(res_unstim, unstim_fileName)

cd3_fileName <- "./Results/diffExpGenes_cd3_all.csv"
write.csv(res_cd3, cd3_fileName)

pma_fileName <- "./Results/diffExpGenes_pma_all.csv"
write.csv(res_pma, pma_fileName)