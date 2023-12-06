################################################################################################################################
################################ RNA Taq-Seq sponge larval condensate experiment ###############################################
################################################################################################################################

library(readxl)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(data.table)
library(vegan)

####---- Input data download: https://apps.aims.gov.au/metadata/view/04c7257a-9941-4378-8fca-2d8dc05e6c67

counts=read.table("allcounts_reordered.txt",header=TRUE,row.names=1)

head(counts) 

length(counts[,1]) #should list the number of genes/isogroups, 13795

names(counts)  #used concentration number in name so coming up as X... e.g. X0A_S94

### assign treatments

Tr <- c(rep("WAF0", 6), rep("WAF2", 6), rep("WAF10", 6), rep("WAF20", 6), rep("WAF50", 6), rep("WAF100", 6))

colData=data.frame(Tr)
head(colData)

### creating a DESeqDataSet

dds<-DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~Tr)

dim(dds) #check to make sure dimensions look correct

### Differential Expression Analysis 

dds<-DESeq(dds)
res<-results(dds)
res

#mcols(res)$description

plotMA(res, main="DESeq2", ylim=c(-2,2))

plotDispEsts(dds)

resOrdered<-res[order(res$padj),] # reorder the results to have the genes with the smallest pval listed first
head(resOrdered)

summary(res) #summary table of results

#write.csv(as.data.frame(resOrdered), file="CartCond_DESeq2results_50vs0.csv")

# cycle through other treatment comparisons - note control/untreated condition needs to be listed last:: log2 fold change (MAP): condition treated vs untreated", meaning that the estimates are of log2(treated / untreated), as would be returned by contrast=c("condition","treated","untreated")
# if your 'control' isn't listed 1st in R alphabetical order you can also set it : dds$condition <- relevel(dds$condition, "untreated")
resTr<-results(dds, contrast=c("Tr", "WAF100","WAF0"))

resTrOrdered<-resTr[order(resTr$padj),]

head(resTrOrdered)

summary(resTr)

#write.csv(as.data.frame(resTrOrdered), file="CartCond_DESeq2results_100vs0.csv")

resTr2<-results(dds, contrast=c("Tr", "WAF20","WAF0"))
resTr2Ordered<-resTr2[order(resTr2$padj),]
head(resTr2Ordered)

summary(resTr2)

#write.csv(as.data.frame(resTr2Ordered), file="CartCond_DESeq2results_20vs0.csv")

resTr3<-results(dds, contrast=c("Tr", "WAF10","WAF0"))
resTr3Ordered<-resTr3[order(resTr3$padj),]
head(resTr3Ordered)

summary(resTr3)

#write.csv(as.data.frame(resTr3Ordered), file="CartCond_DESeq2results_10vs0.csv")

resTr4<-results(dds, contrast=c("Tr", "WAF2","WAF0"))
resTr4Ordered<-resTr4[order(resTr4$padj),]
head(resTr4Ordered)

summary(resTr4)

#write.csv(as.data.frame(resTr4Ordered), file="CartCond_DESeq2results_2vs0.csv")

plotCounts(dds, gene=which.min(res$padj), intgroup="Tr") #plots the most significant isogroup14113_c0_g1, normalized counts (WAF50 vs 0)
plotCounts(dds, gene=which.min(resTr$padj), intgroup="Tr") #plots the most significant isogroup7641_c0_g1, normalized counts (WAF100 vs 0)
plotCounts(dds, gene=which.min(resTr2$padj), intgroup="Tr") # isogroup38837_c0_g2, normalized counts (WAF20 vs 0)
plotCounts(dds, gene=which.min(resTr3$padj), intgroup="Tr") #isogroup38837_c0_g2, normalized counts (WAF10 vs 2)
plotCounts(dds, gene=which.min(resTr4$padj), intgroup="Tr") #isogroup4771_c0_g2, normalized counts (WAF2 vs 0)

# Make a column with the isogroup labels from the rownames

DF.res <- as.data.frame(cbind("isogroup"=row.names(res), res))
DF.resTr <- as.data.frame(cbind("isogroup"=row.names(resTr), resTr))
DF.resTr2 <- as.data.frame(cbind("isogroup"=row.names(resTr2), resTr2))
DF.resTr3 <- as.data.frame(cbind("isogroup"=row.names(resTr3), resTr3))
DF.resTr4 <- as.data.frame(cbind("isogroup"=row.names(resTr4), resTr4))

# Rename the columns to prevent merging of columns that vary among the different comparisons
# Only isogroup and baseMean columns are the same

colnames(DF.res) <- c("isogroup", "baseMean", "res_log2FoldChange", "res_lfcSE", "res_stat", "res_pvalue", "res_padj")
colnames(DF.resTr) <- c("isogroup", "baseMean", "resTr_log2FoldChange", "resTr_lfcSE", "resTr_stat", "resTr_pvalue", "resTr_padj")
colnames(DF.resTr2) <- c("isogroup", "baseMean", "resTr2_log2FoldChange", "resTr2_lfcSE", "resTr2_stat", "resTr2_pvalue", "resTr2_padj")
colnames(DF.resTr3) <- c("isogroup", "baseMean", "resTr3_log2FoldChange", "resTr3_lfcSE", "resTr3_stat", "resTr3_pvalue", "resTr3_padj")
colnames(DF.resTr4) <- c("isogroup", "baseMean", "resTr4_log2FoldChange", "resTr4_lfcSE", "resTr4_stat", "resTr4_pvalue", "resTr4_padj")

list.results <- list(DF.res, DF.resTr, DF.resTr2, DF.resTr3, DF.resTr4)
all.treatment.adjusted.pvalues <- Reduce(function(...) merge(..., all=T), list.results)
write.csv(as.data.frame(all.treatment.adjusted.pvalues), file="DESeq2results_combined.csv")

# Append gene names to full results dataset
colnames(CartGene) <- c("isogroup", "Name")
list.results2 <- list(CartGene, all.treatment.adjusted.pvalues)
all.treatment.adjusted.pvalues_wGeneName <- Reduce(function(...) merge(..., all=T), list.results2)
write.csv(as.data.frame(all.treatment.adjusted.pvalues_wGeneName), file="DESeq2results_combined_wGeneNames.csv")

## Venn diagram

#p=read.csv("CartCond_combined_pvals.csv", row.names = 1) #dont forget row.names otherwise it will pull out row numbers and not gene names below

p <- all.treatment.adjusted.pvalues[,c(1,7,12,17,22, 27)]
row.names(p) <- p[,1]
colnames(p) <- c("isogroup", "padj.50", "padj.100", "padj.20", "padj.10", "padj.2")

# padj 0.05 cutoff

WAF2=row.names(p[p$padj.2<=0.05 & !is.na(p$padj.2),])
WAF10=row.names(p[p$padj.10<=0.05 & !is.na(p$padj.10),])
WAF20=row.names(p[p$padj.20<=0.05 & !is.na(p$padj.20),])
WAF50=row.names(p[p$padj.50<=0.05 & !is.na(p$padj.50),])
WAF100=row.names(p[p$padj.100<=0.05 & !is.na(p$padj.100),])

candidates=list("2% WAF" = WAF2, "10% WAF" = WAF10, "20% WAF" = WAF20, "50% WAF" = WAF50, "100% WAF" = WAF100)

venn(candidates)   

tmp <- venn(candidates, show.plot = FALSE)

tmp.inter <- attr(tmp, "intersections")
tmp.inter

#capture.output(tmp.inter, file = "IsogroupsIDs_Venn.txt") #still would need editing in excel

# alternative - each vennID corresponds to headings in tmp.inter    

Venn1 <- as.data.frame(cbind("Isogroup" = tmp.inter$`2% WAF:10% WAF`, "venn_1_ids" = 1))
Venn2 <- as.data.frame(cbind("Isogroup" = tmp.inter$`2% WAF:100% WAF`, "venn_2_ids" = 1))
Venn3 <- as.data.frame(cbind("Isogroup" = tmp.inter$`20% WAF:50% WAF`, "venn_3_ids" = 1))
Venn4 <- as.data.frame(cbind("Isogroup" = tmp.inter$`20% WAF:100% WAF`, "venn_4_ids" = 1))
Venn5 <- as.data.frame(cbind("Isogroup" = tmp.inter$`50% WAF:100% WAF`, "venn_5_ids" = 1))
Venn6 <- as.data.frame(cbind("Isogroup" = tmp.inter$`2% WAF:20% WAF:100% WAF`, "venn_6_ids" = 1))
Venn7 <- as.data.frame(cbind("Isogroup" = tmp.inter$`20% WAF:50% WAF:100% WAF`, "venn_7_ids" = 1))
Venn8 <- as.data.frame(cbind("Isogroup" = tmp.inter$`10% WAF:20% WAF:50% WAF`, "venn_8_ids" = 1))
Venn9 <- as.data.frame(cbind("Isogroup" = tmp.inter$`2% WAF`, "venn_9_ids" = 1))
Venn10 <- as.data.frame(cbind("Isogroup" = tmp.inter$`10% WAF`, "venn_10_ids" = 1))
Venn11 <- as.data.frame(cbind("Isogroup" = tmp.inter$`20% WAF`, "venn_11_ids" = 1))
Venn12 <- as.data.frame(cbind("Isogroup" = tmp.inter$`50% WAF`, "venn_12_ids" = 1))
Venn13 <- as.data.frame(cbind("Isogroup" = tmp.inter$`100% WAF`, "venn_13_ids" = 1))

list.venn <- list(Venn1, Venn2, Venn3, Venn4, Venn5, Venn6, Venn7, Venn8, Venn9, Venn10, Venn11, Venn12, Venn13) 
#capture.output(list.venn, file = "VennList.txt")
all.venn.ids <- Reduce(function(...) merge(..., all=T), list.venn)
write.table(all.venn.ids, file="IsogroupIDs_Venn.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#------log2FC data for genes of interest: stress response network & Nrf2------

p2 <- all.treatment.adjusted.pvalues[,c(1,3,8,13,18,23)]
row.names(p2) <- p2[,1]
colnames(p2) <- c("isogroup", "log2FC.50", "log2FC.100", "log2FC.20", "log2FC.10", "log2FC.2")
# reorder columns from 2 to 100 WAF
p2 <- p2[,c(1,6,5,4,2,3)]

goi <- read.csv(file = "Genes_of_interest.csv", header = TRUE)

goi.log2FC <- inner_join(goi, p2, by = "isogroup")

write.csv(goi.log2FC, file = "log2FC_genes_of_interest.csv")

#------add p-values and isogroup names to above

goi.log2FC_pval <- inner_join(goi.log2FC, p, by = 'isogroup')

colnames(goi.log2FC_pval) <- c("Isogroup",   "Gene.ID",   "Gene.Name",  "Gene.group", "log2FC.2",   "log2FC.10",  "log2FC.20",  "log2FC.50",  "log2FC.100", "padj.50", "padj.100",   "padj.20", "padj.10",  "padj.2") 

goi.log2FC_pvalName <- inner_join(goi.log2FC_pval, CartGene, by = "Isogroup")

goi.log2FC_pvalName <- goi.log2FC_pvalName %>% select(Isogroup, Gene.ID, Gene.Name, Gene.group, Name, log2FC.2, padj.2, log2FC.10, padj.10, log2FC.20, padj.20, log2FC.50, padj.50, log2FC.100, padj.100)

write.csv(goi.log2FC_pvalName, file = "genes.of.interest_AllDESeq2dat.csv")

#--------log2FC data for genes of interest: GO:0006979 response to oxidative stress

goi2 <- read.csv(file = "GO_0006979_genes.csv", header = TRUE)

goi2.log2FC <- inner_join(goi2, p2, by = "isogroup")

colnames(goi2.log2FC) <- c("Isogroup", "log2FC.2", "log2FC.10", "log2FC.20", "log2FC.50", "log2FC.100")

goi2.log2FC.names <- inner_join(goi2.log2FC, CartGene, by = "Isogroup")

write.csv(goi2.log2FC.names, file = "log2FC_GO_0006979_wGeneNames.csv")

##########################################################
####-----Data transformations and visualization -----#####
##########################################################

# varianceStabilizingTransformation

vsd <- varianceStabilizingTransformation(dds)

vstMat <- assay(vsd)

#write.csv(vstMat, file = "CartCond_VSTdata.csv")


#select the top 30 genes to plot in heatmap

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))

heatmap.2(assay(rlog)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = TRUE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

# rescale vsd 'select' expression data

meanVSD <- apply(assay(vsd)[select,], 1, mean)  # select: 30, 50,  100
meanVSDmat <- assay(vsd)[select,]-meanVSD

hmcol2 <- colorRampPalette(brewer.pal(11, "PuOr"))(100)

heatmap.2(meanVSDmat, col = hmcol2,
          Rowv = TRUE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

heatmap.2(meanVSDmat, col=hmcol2, Rowv = TRUE, Colv = FALSE, dendrogram = c("none"), scale="none", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA, 
          ColSideColors = c(rep("#2b4b9b", 6), rep("#8ebe21",6), rep("#a6529a", 6), rep("#66c1bf", 6), rep("#ef7e1b", 6), rep("#483c32", 6)))
legend("bottomleft",
       legend=c("0% WAF", "2% WAF", "10% WAF", "20% WAF", "50% WAF", "100% WAF"), # category labels
       fill=c("#2b4b9b", "#8ebe21", "#a6529a", "#66c1bf", "#ef7e1b", "#483c32"),
       border=FALSE,          # doesn't plot lines around the treatment colors in the legend
       bty="n"    # doesn't plot a box around the legend
       # y.intersp = 0.7,          # character interspacing factor for vertical (y) line distances
       # cex=0.7          # makes the font smaller
)

distsRL <- dist(t(assay(rlog)))
mat <- as.matrix(distsRL)
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))

# Convert the DESeq / vst transformed objects to a data frames see: https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/

# 100 vs 0 = resTr

resTr.dat <- as.data.frame(resTr)

head(resTr.dat)

# count number of significant genes 

sum(resTr.dat$padj < 0.05, na.rm = TRUE)  

##############################
### bar plot of up and down regulated genes for all treatments ###

DE.barplot.dat <- read_excel("Combined DEseq2 data.xlsx", sheet = "NoDEsUpDown")

ggplot(DE.barplot.dat, aes(Tr), ylim(-515:420)) + 
  geom_bar(data = subset(DE.barplot.dat, variable == "count.up"), 
           aes(y = value, fill = DE), stat = "identity", position = "dodge") +
  geom_bar(data = subset(DE.barplot.dat, variable == "count.down"), 
           aes(y = -value, fill = DE), stat = "identity", position = "dodge") + 
  geom_hline(yintercept = 0,colour = "grey90") +
  scale_fill_manual(values = c("#F3A446", "#9F97C4"))+
  theme_classic()


last_plot() + 
  geom_text(data = subset(DE.barplot.dat, variable == "count.up"), 
            aes(Tr, value, group=DE, label=value),
            position = position_dodge(width=0.9), vjust = 1.5, size=4) +
  geom_text(data = subset(DE.barplot.dat, variable == "count.down"), 
            aes(Tr, -value, group=DE, label=value),
            position = position_dodge(width=0.9), vjust = -.5, size=4) +
  coord_cartesian(ylim = c(-515, 420))

### need to remove NA from padj to get below code to work..


resTr.dat$sig <- ifelse(resTr.dat$padj < 0.1, "sig", na.rm = TRUE)


vstMat2 <- as.data.frame(vstMat)
vstMat2$Isogroup <- rownames(vstMat)
head(vstMat2)


### subset significantly expressed genes for 100 vs 0 with log2fc > 1.5

sig100 <- rownames(resTr.dat[resTr.dat$padj <= .05 & abs(resTr.dat$log2FoldChange) > 1.5,])

head(sig100)

vstMat100 <- vstMat2[vstMat2$Isogroup %in% sig100,]

###### convert to long format for ggplot2 heatmap ****** I like heatmap.2 better so need to convert vstMat100 back to matrix ******

##vstMat100_long <- reshape2::melt(vstMat100, id.vars=c("Isogroup"))

##head(vstMat100_long)

##library(viridis)

##ggplot(vstMat100_long, aes(x=variable, y=Isogroup, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())

vstMat100.mat <- select(vstMat100, c(1:36)) #get rid of Isogroup column added before

# turn back into a matrix

vstMat100.mat <- data.matrix(vstMat100.mat)

heatmap.2(vstMat100.mat, col=hmcol2, Rowv = TRUE, Colv = FALSE, dendrogram = c("row"), scale="none", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)

# rescale up and down
meanVSTmat100 <- apply(vstMat100.mat, 1, mean) 
vstMat100.mat2 <- vstMat100.mat-meanVSTmat100

###### Final publication heatmap ##### 
heatmap.2(vstMat100.mat2, col=hmcol2, Rowv = TRUE, Colv = FALSE, dendrogram = c("none"), scale="none", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA, 
          ColSideColors = c(rep("#2b4b9b", 6), rep("#8ebe21",6), rep("#a6529a", 6), rep("#66c1bf", 6), rep("#ef7e1b", 6), rep("#483c32", 6)))
legend("bottomleft",
       legend=c("0% WAF", "2% WAF", "10% WAF", "20% WAF", "50% WAF", "100% WAF"), # category labels
       fill=c("#2b4b9b", "#8ebe21", "#a6529a", "#66c1bf", "#ef7e1b", "#483c32"),
       border=FALSE,          # doesn't plot lines around the treatment colors in the legend
       bty="n"    # doesn't plot a box around the legend
       # y.intersp = 0.7,          # character interspacing factor for vertical (y) line distances
       # cex=0.7          # makes the font smaller
)

####---- PCA ----####

plotPCA(rlog, intgroup=c("Tr"))

pca.dat <- plotPCA(vsd, intgroup=c("Tr"), returnData = TRUE)

percentVar <- round(100 * attr(pca.dat, "percentVar"))

pca.plot <- ggplot(pca.dat, aes(PC1, PC2, color=Tr)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

####---- Final publication PCA ----####

# reorder treatments on plot
newtab = data.table(pca.plot$data)
newtab$Tr <- ordered(newtab$Tr, levels = c("WAF0","WAF2","WAF10", "WAF20", "WAF50", "WAF100"))
pca.plot$data <- newtab
print(pca.plot)
#pca.plot + scale_color_hue()
plot.colors <- c("#0000ff","#80ff00", "#ff00ff", "#40e0d0", "#ff8000", "#483c32")
pca.plot + scale_color_manual(values=plot.colors) + theme_classic()

######################## Dispersion test & PERMANOVA ########################

metadata <- read.csv(file = "metadata_permanova.csv", header = TRUE, sep = ",")

metadata$treatment <- as.factor(metadata$treatment)

vstMat.t <- t(vstMat[,1:36])

dist <- vegdist(vstMat.t, method = "bray")

beta <- betadisper(dist, metadata$treatment)

disper.test <-  permutest(beta, permutations =9999) # pval = 0.48 i.e., > 0.05 so OK

perm <- adonis2(dist ~ treatment,
        data = metadata, permutations = 9999)

perm

pairwise <- pairwise.adonis2(dist ~ treatment, data = metadata, p.adjust.methods = "fdr")

pairwise 

