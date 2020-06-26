setwd("~/GitHub/stroke-trf")
rm(list = ls())

pacman::p_load(DESeq2, apeglm, topGO)

options(stringsAsFactors = F)

get.count.change <- function(BM, LF) {
  return((BM * 2^LF) - BM)
}

remove.char <- function(str, n) {
  return(substr(str, n+1, nchar(str)))
}

# DIFFERENTIAL EXPRESSION OF LARGE RNAs ####
#load data
lnames<- load(file = "data/wgcna_input_gene.RData");
lnames

#compile metadata
sample.table <- traitData[, c("disease.state", "age", "batch", "sample_id", "mir_filename")]
meta <- read.csv("data/raw/meta.csv", sep = ";", dec = ",")

sample.table$time2blood <- meta$time2blood[match(sample.table$mir_filename, meta$FileName)]


#factorization binning of continuous variables
sample.table$age <- cut(sample.table$age, 5)
sample.table$time2blood <- cut(sample.table$time2blood, 5)

set.seed(1)
sample.table$time2blood[is.na(sample.table$time2blood)] <- sample(sample.table$time2blood[!is.na(sample.table$time2blood)],
                                                                  length(sample.table$time2blood[is.na(sample.table$time2blood)]))


#gene DESeq####
counts <- read.csv("data/raw/mrna_salmon_counts.csv", row.names = 1, check.names = F)
head(counts)

## >sample table ####
dim(counts)
dim(sample.table)
colnames(counts) == rownames(sample.table)
all(colnames(counts) %in% rownames(sample.table))
counts <- counts[, rownames(sample.table)]

## >DESeq ####
se <- SummarizedExperiment(as.matrix(round(counts, 0)), 
                           colData = sample.table)
ddsSE <- DESeqDataSet(se, design = ~ age + batch + disease.state)
ddsSE

#>DESeq run ####
dds <- estimateSizeFactors(ddsSE)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=500)
resultsNames(dds)
#stroke vs control
stroke <- lfcShrink(dds, coef = c("disease.state_stroke_vs_control"), type = "apeglm")
plotMA(stroke)
summary(stroke)

ensg2symbol <- readRDS(file = "data/ensg2symbol_all.rds")
stroke$symbol <- ensg2symbol$hgnc_symbol[match(rownames(stroke), ensg2symbol$ensembl_gene_id)]
stroke$countChange <- get.count.change(stroke$baseMean, stroke$log2FoldChange)
stroke$ensg <- rownames(stroke)
stroke <- stroke[order(stroke$padj),]
stroke <- stroke[!is.na(stroke$padj),]
stroke

write.csv(stroke, file = "out/DESeq2_gene_stroke_vs_con.csv",
          row.names = F, quote = F)

strokeSig <- subset(stroke, padj < .05)
#names
strokeSig$ensg <- rownames(strokeSig)
strokeSig$countChange <- get.count.change(strokeSig$baseMean, strokeSig$log2FoldChange)
(strokeSig <- strokeSig[order(strokeSig$padj, decreasing = F),])
(strokeSig <- strokeSig[order(abs(strokeSig$countChange), decreasing = T),])
as.data.frame(strokeSig)
nrow(strokeSig) #number DE

nrow(strokeSig[strokeSig$log2FoldChange>0,])
nrow(strokeSig[strokeSig$log2FoldChange<0,])
nrow(strokeSig[strokeSig$log2FoldChange>1.4,])
nrow(strokeSig[strokeSig$log2FoldChange<(-1.4),])

#plot####
dev.off()
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

hist(stroke$log2FoldChange, xlim = c(-1,1), breaks = 1000)

topT <- as.data.frame(stroke)

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", #col = color,
                cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<0.05 & abs(log2FoldChange)>1.4), text(log2FoldChange, -log10(padj), labels=subset(topT$symbol, topT$padj<0.05 & abs(topT$log2FoldChange)>1.4), cex=0.8, pos=3))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1.4, col="black", lty=4, lwd=2.0)
abline(v=1.4, col="black", lty=4, lwd=2.0)
abline(h=-log10(.05), col="black", lty=4, lwd=2.0)

#GO ANALYSIS####
geneID2GO <- readRDS(file = "data/geneID2GO_genes.rds")
#>order by pval####
stroke <- stroke[!is.na(stroke$padj),]

#>>all####
genes <- stroke[order(stroke$padj),]
#ontology
ontology <- "BP"
#description
description <- "top genes DE in stroke"
#allGenes list
allGenes <- genes$padj[1:2000]
names(allGenes) <- genes$ensg[1:2000]
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore < .05)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   weight = resultWeight, 
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
(sig <- allRes[allRes$weight<.05,])
write.csv(sig, file = "out/gene_GO_significant_all.csv")

#>>>log fold change > 1.4####
genes <- stroke[order(stroke$padj),]
idx <- which(abs(stroke$log2FoldChange) < 1.4 & stroke$padj < .05)
genes <- genes[-idx,]
#ontology
ontology <- "BP"
#description
description <- "top genes DE in stroke"
#allGenes list
allGenes <- genes$padj[1:2000]
names(allGenes) <- genes$ensg[1:2000]
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore < .05)
}

#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   weight = resultWeight, 
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
(sig <- allRes[allRes$weight<.05,])
write.csv(sig, file = "out/gene_GO_significant_all_above_1.4.csv")


#>>positive####
genes <- stroke[stroke$log2FoldChange>0,]
genes <- genes[order(genes$padj),]
#ontology
ontology <- "BP"
#description
description <- "top genes DE in stroke"
#allGenes list
allGenes <- genes$padj[1:2000]
names(allGenes) <- genes$ensg[1:2000]
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore < .05)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   weight = resultWeight,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
(sig <- allRes[allRes$weight<.05,])
write.csv(sig, file = "out/gene_GO_significant_pos.csv")

#>>>1.4####
genes <- stroke[stroke$log2FoldChange>0,]
genes <- genes[order(genes$padj),]
idx <- which(abs(stroke$log2FoldChange) < 1.4 & stroke$padj < .05)
genes <- genes[-idx,]
#ontology
ontology <- "BP"
#description
description <- "top genes DE in stroke"
#allGenes list
allGenes <- genes$padj[1:2000]
names(allGenes) <- genes$ensg[1:2000]
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore < .05)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   weight = resultWeight,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
(sig <- allRes[allRes$weight<.05,])
write.csv(sig, file = "out/gene_GO_significant_pos_above_1.4.csv")

#>>negative####
genes <- stroke[stroke$log2FoldChange<0,]
genes <- genes[order(genes$padj),]
#ontology
ontology <- "BP"
#description
description <- "top genes DE in stroke"
#allGenes list
allGenes <- genes$padj[1:2000]
names(allGenes) <- genes$ensg[1:2000]
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore < .05)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   weight = resultWeight,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
(sig <- allRes[allRes$weight<.05,])
write.csv(sig, file = "out/gene_GO_significant_neg.csv")

#>>>1.4####
genes <- stroke[stroke$log2FoldChange<0,]
genes <- genes[order(genes$padj),]
idx <- which(abs(stroke$log2FoldChange) < 1.4 & stroke$padj < .05)
genes <- genes[-idx,]
#ontology
ontology <- "BP"
#description
description <- "top genes DE in stroke"
#allGenes list
allGenes <- genes$padj[1:2000]
names(allGenes) <- genes$ensg[1:2000]
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore < .05)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   weight = resultWeight,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
(sig <- allRes[allRes$weight<.05,])
write.csv(sig, file = "out/gene_GO_significant_neg_above_1.4.csv")


#>countChange####
stroke <- stroke[!is.na(stroke$countChange),]
mean(genes$ensg %in% names(geneID2GO))

#>>positive####
genes <- stroke[order(stroke$padj < .05, stroke$countChange, decreasing = T),]
num_de <- sum(genes$padj<.05 & genes$countChange > 0)
top <- 100
#ontology
ontology <- "BP"
#description
description <- "top genes DE in stroke"
#allGenes list
allGenes <- genes$countChange[1:2000]
names(allGenes) <- genes$ensg[1:2000]
allGenes[(top+1):2000] <- 1
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore > 1)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   weight = resultWeight,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
(sig <- allRes[allRes$weight<.05,])
write.csv(sig, file = "out/gene_GO_top_100_pos.csv")


#>>negative####
genes <- stroke[order(stroke$padj < .05, stroke$countChange, decreasing = F),]
num_de <- sum(genes$padj<.05 & genes$countChange < 0)
top <- 100
#ontology
ontology <- "BP"
#description
description <- "top genes DE in stroke"
#allGenes list
allGenes <- genes$countChange[1:2000]
names(allGenes) <- genes$ensg[1:2000]
allGenes[(top+1):2000] <- 1
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore < 1)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   weight = resultWeight,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
(sig <- allRes[allRes$weight<.05,])
write.csv(sig, file = "out/gene_GO_top_100_neg.csv")


#>>all####
genes <- stroke[order(stroke$padj < .05, abs(stroke$countChange), decreasing = T),]
num_de <- sum(genes$padj<.05)
top <- 100
#ontology
ontology <- "BP"
#description
description <- "top genes DE in stroke"
#allGenes list
allGenes <- abs(genes$countChange)[1:2000]
names(allGenes) <- genes$ensg[1:2000]
allGenes[(top+1):2000] <- 1
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore > 1)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   weight = resultWeight,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
(sig <- allRes[allRes$weight<.05,])
write.csv(sig, file = "out/gene_GO_top_100_abs.csv")
