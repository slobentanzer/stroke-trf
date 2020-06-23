setwd("~/GitHub/stroke-trf")
rm(list = ls())

library(DESeq2)
library(apeglm)
library(ggplot2)
library(waffle)
library(RColorBrewer)
library(Biostrings)
library(BiocParallel)
library(Rtsne)
library(pheatmap)

get.count.change <- function(BM, LF) {
  return((BM * 2^LF) - BM)
}

remove.char <- function(str, n) {
  return(substr(str, n+1, nchar(str)))
}

# DIFFERENTIAL EXPRESSION OF SMALL RNAs ####

#load data
tRFInfo <- readRDS("data/tRFInfo_human.rds")
lnames <- load(file = "data/cholinergic_associated_smrnas.RData")
lnames

lnames<- load(file = "data/wgcna_input_smrna.RData");
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

#Read counts and unify
counts_t <- read.csv("data/raw/tRNA_reads_stroke.csv", row.names = 1)
colnames(counts_t) <- gsub("X", "", colnames(counts_t))
counts_t[is.na(counts_t)] <- 0
dim(counts_t);
names(counts_t);

counts_m <- read.csv("data/raw/miRNA_reads_stroke.csv", row.names = 1)
colnames(counts_m) <- gsub("X", "", colnames(counts_m))
counts_m[is.na(counts_m)] <- 0
dim(counts_m);
names(counts_m);

colnames(counts_m) == colnames(counts_t)
colnames(counts_m) %in% colnames(counts_t)

counts_t <- counts_t[, rownames(sample.table)]
counts_m <- counts_m[, sample.table$mir_filename]

rownames(sample.table) == colnames(counts_t)
colnames(counts_m) <- sample.table$sample_id #use sample ids
colnames(counts_t) <- sample.table$sample_id
rownames(sample.table) <- sample.table$sample_id

counts <- rbind(counts_m, counts_t)

#smRNA DESeq####
## >DESeq ####
se <- SummarizedExperiment(as.matrix(counts), 
                           colData = sample.table)
ddsSE <- DESeqDataSet(se, design = ~ age + batch + time2blood + disease.state)
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
stroke$name <- rownames(stroke)
stroke$countChange <- get.count.change(stroke$baseMean, stroke$log2FoldChange)
stroke <- stroke[order(abs(stroke$countChange), decreasing = T),]
stroke <- stroke[!is.na(stroke$padj),]

write.csv(stroke, file = "out/DESeq2_smrna_stroke_vs_con.csv", row.names = F, quote = F)

strokeSig <- data.frame(subset(stroke, padj < .05))

#amounts
strokeSig_mir <- strokeSig[-grep("tRF", strokeSig$name),]
strokeSig_trf <- strokeSig[grep("tRF", strokeSig$name),]

#tRF information####
strokeSig_trf$aa <- tRFInfo$`Amino Acid`[match(strokeSig_trf$name, tRFInfo$MINTplate)]
strokeSig_trf$type <- tRFInfo$`tRF type(s)`[match(strokeSig_trf$name, tRFInfo$MINTplate)]

table(strokeSig_trf$type)

table(strokeSig_trf$aa)

parts <- table(strokeSig_trf$aa)
parts <- parts[order(parts, decreasing = T)]

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
qual_col_pals <- qual_col_pals[c("Dark2", "Set1", "Accent", "Set2", "Set3", "Pastel1", "Pastel2", "Paired"),]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
waffle(parts, colors = col_vector[1:length(parts)], rows = 12, size = 1)
ggsave("img/tRF_amino_acid_waffle_diagram.svg", width = 8, height = 5)

#plot volcanoes####
####trf
dev.off()
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

topT <- as.data.frame(stroke[grep("tRF", rownames(stroke)),])
topT$color <- "darkgrey"

topT$color[topT$name %in% cholinotrfs] <- "darkred"
sort(table(topT$color))

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col = alpha(color, .7),
                cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1.4, col="black", lty=4, lwd=2.0)
abline(v=1.4, col="black", lty=4, lwd=2.0)
abline(h=-log10(.05), col="black", lty=4, lwd=2.0)

svg("img/trf_volcano_mod_color.svg", width = 8, height = 8)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col = alpha(color, .7),
                cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1.4, col="black", lty=4, lwd=2.0)
abline(v=1.4, col="black", lty=4, lwd=2.0)
abline(h=-log10(.05), col="black", lty=4, lwd=2.0)
dev.off()


####mir
dev.off()
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

topT <- as.data.frame(stroke[-grep("tRF", rownames(stroke)),])
topT$color <- "darkgrey"

topT$color[topT$name %in% cholinomirs] <- "darkred"
sort(table(topT$color))

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col = alpha(color, .7),
                cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1.4, col="black", lty=4, lwd=2.0)
abline(v=1.4, col="black", lty=4, lwd=2.0)
abline(h=-log10(.05), col="black", lty=4, lwd=2.0)

svg("img/mir_volcano_mod_color.svg", width = 8, height = 8)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col = alpha(color, .7),
                cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1.4, col="black", lty=4, lwd=2.0)
abline(v=1.4, col="black", lty=4, lwd=2.0)
abline(h=-log10(.05), col="black", lty=4, lwd=2.0)
dev.off()

#save####
save(stroke, 
  strokeSig, 
  file = "data/DE_sig.RData")
load("data/DE_sig.RData")


#HOMOLOGY FOR ALL DE#####
#compute for tRFs in variability filter
tRFs <- rownames(strokeSig) #use regular DE without interaction terms
tRFs <- tRFs[grep("tRF", tRFs)]
write.table(tRFs, "out/de_trfs.txt", sep = "\n", col.names = F, row.names = F, quote = F)

### --- CONVERT PLATES TO SEQUENCES USING MINTPLATE ALGORITHM --- ###

tRFs <- read.table("data/de_seqs.txt", fill = T, skip = 5)
tRFs <- tRFs[1:(grep("Thank", tRFs$V1)-1),]
colnames(tRFs) <- c("seq", "MINTplate")

get.homology <- function(nam){
  if(!file.exists(paste0("data/temp/", nam, ".rds"))){
    seq <- tRFs$seq[tRFs$MINTplate == nam]
    lst <- vector(mode = "numeric", length = nrow(tRFs))
    for(i in 1:nrow(tRFs)) {
      # if(i%%1000 == 0)
      #   message(i)
      lst[i] <- pairwiseAlignment(seq, tRFs$seq[i], type = "global", scoreOnly = T)
    }
    saveRDS(lst, file = paste0("data/temp/", nam, ".rds"))
  }
}

# use DE tRFs for input
run <- F #only run once, takes time
if(run){
  TOI <- rownames(strokeSig[grepl("tRF", rownames(strokeSig)), ])
  fils <- list.files(paste0("data/temp/"))
  fils <- gsub(".rds", "", fixed = T, fils)
  TOI <- TOI[!TOI %in% fils]
  bplapply(TOI, get.homology, BPPARAM = MulticoreParam(workers = 7))
  #read rds
  dir <- "data/temp/"
  files <- paste0(dir, list.files(dir))
  homologs <- lapply(files, function(file) {
    return(readRDS(file))
  })
  names(homologs) <-
    gsub("data/temp/", "", gsub(".rds", "", fixed = T, files))
  homologs <- plyr::ldply(homologs, id = T)
  rownames(homologs) <- homologs$.id
  homologs$.id <- NULL
  colnames(homologs) <- tRFs$MINTplate
  homologs[1:5, 1:5]
  dim(homologs)
  saveRDS(homologs, file = "data/trf_homologs_hsa.rds")
}
#heatmap####
homologs <- readRDS(file = "data/trf_homologs_hsa.rds")
pheatmap(homologs)

#binary HM for categorisation
fam <- data.frame(apply(homologs, 2, function(x) {
  vec <- x
  vec[which(vec<=0)] <- 0
  vec[which(vec>0)] <- 1 #optional for more clarity
  return(as.numeric(vec))
}))
rownames(fam) <- rownames(homologs)

heat <- pheatmap(as.matrix(fam),
                 # color = c("grey60", "steelblue"),
                 border_color = "black"#,
                 # cellwidth = 10, cellheight = 10,
                 # filename = "img/HR0V_SPOX_homology_binary_hsa.png"
)
rows <- heat$tree_row$order

#continuous HM for plotting
fam <- data.frame(apply(homologs, 2, function(x) {
  vec <- x
  # vec[which(vec<=0)] <- 0
  # vec[which(vec>0)] <- 1 #optional for more clarity
  return(vec)
}))
rownames(fam) <- rownames(homologs)

fam <- fam[rows, ]

# make color palette centered at 0####
# Returns a vector of 'num.colors.in.palette'+1 colors. The first 'cutoff.fraction'
# fraction of the palette interpolates between colors[1] and colors[2], the remainder
# between colors[3] and colors[4]. 'num.colors.in.palette' must be sufficiently large
# to get smooth color gradients.
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) %% 2 == 0)
  ramp1 <- colorRampPalette(colors[1:(length(colors)/2)], interpolate = "linear")(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[(length(colors)/2+1):length(colors)], interpolate = "linear")(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}

dist <- seq(min(fam), max(fam), length = 100)
cutoff.distance <- ecdf(dist)(0)
blu <- c(RColorBrewer::brewer.pal(9, "Blues"))
red <- c(RColorBrewer::brewer.pal(9, "YlOrRd"))
cols <- makeColorRampPalette(c(rev(blu), red),
                             cutoff.distance,
                             100)
#plot large HM
pheatmap(as.matrix(fam),
         color = cols,
         border_color = "black",
         cluster_rows = F,
         cellwidth = 10, cellheight = 10,
         filename = "img/DE_tRF_homology_hsa.pdf"
); dev.off()

#homologies greater zero
idx <- which(apply(fam, 2, function(x) any(x>0)))
fam_gz <- fam[, idx]
#plot
p <- pheatmap(as.matrix(fam_gz),
         color = cols,
         border_color = "black"#,
         # cluster_rows = F,
         # cellwidth = 10, cellheight = 10,
         # filename = "img/DE_tRF_homology_hsa_only_greater_zero.pdf"
)
dim(fam_gz)

#>display via tSNE####

cols <- col_vector

#zoom, text####
i <- 16
set.seed(1234)
cl <- cutree(p$tree_row, i)
tsne <- Rtsne(as.matrix(fam_gz))
svg(paste0("img/tsne", i, "_zoom_label.svg"), width = 10, height = 10)
plot(tsne$Y, col = cols[cl], pch = 16, main = i)
text(tsne$Y[,1], y = tsne$Y[,2], labels = names(cl))
dev.off()

#color by amino acid origin####
set.seed(1234)
tsne <- Rtsne(as.matrix(fam_gz),
              perplexity = i,
              max_iter = 5000)
names <- rownames(fam_gz)
aa <-
  factor(tRFInfo$`Amino Acid`[match(names, tRFInfo$MINTplate)], levels = names(parts)) #from waffle diagram, to get same colours!
svg(
  paste0("img/tsne_DE_tRFs_by_amino_acid_", i, ".svg"),
  width = 6,
  height = 6
)
plot(tsne$Y,
     col = cols[aa],
     # xlim = c(-11, 11),
     pch = 16,
     main = paste0("perp: ", i))
# text(tsne$Y[,1] + .2, y = tsne$Y[,2] - .1, labels = aa)
legend(
  -20,
  15,
  horiz = F,
  ncol = 4,
  legend = names(parts),
  col = cols,
  pch = 16
)
dev.off()

