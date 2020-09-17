setwd("~/GitHub/stroke-trf")
rm(list = ls())

pacman::p_load(
  tidyverse,
  pheatmap,
  fitdistrplus,
  DESeq2,
  ggplot2,
  Rtsne,
  RColorBrewer,
  ggrepel
  )

options(stringsAsFactors = FALSE);

#Read in
meta <- read.csv("data/raw/juzenas_sample_meta.csv", sep = ",", row.names = 1)
meta$cells <- factor(meta$cells)
head(meta)
rownames(meta) <- meta$file
meta$batch <- factor(meta$batch)

counts <- read.csv("data/raw/juzenas_trna_exclusive_combined_data.csv", sep = ",", row.names = 1)
head(counts)
nrow(counts)
mean(rowMeans(counts) == 0)
mean(colMeans(counts) == 0)

counts <- counts[,colMeans(counts) != 0]

meta <- meta[match(colnames(counts), meta$file),]

#deseq normalize VST
se <- SummarizedExperiment(as.matrix(counts), 
                           colData = meta)
dds <- DESeqDataSet(se, design = ~ cells)

cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)

#FORMALIZE THE DEFINITION OF PRESENCE/ABSENCE####
#collect counts per cell type in a list
if(!file.exists("data/counts_per_cell_list.rds")){
  dis <- vector(mode = "list", length = nrow(cts))
  names(dis) <- rownames(cts)
  for (trf in rownames(counts)) {
    row <- counts[trf,]
    expr <- logical()
    lis <- vector(mode = "list", length = length(unique(meta$cells)))
    names(lis) <- unique(meta$cells)
    for (cell in unique(meta$cells)) {
      lis[[cell]] <- as.numeric(row[meta$file[meta$cells == cell]])
    }
    dis[[trf]] <- lis
  }
  saveRDS(dis, "data/counts_per_cell_list.rds")
} else {
  dis <- readRDS("data/counts_per_cell_list.rds")
}

#>implement statistical test####
cell_types <- levels(meta$cells)

#kolmogorov-smirnov over log-normal distribution of counts with estimated mean and sd
if(!file.exists("data/pvals_per_cell_ks_list.rds")){
  pvals <- vector(mode = "list", length = length(dis))
  names(pvals) <- names(dis)
  
  
  for(trf in names(dis)) {
    p.temp <- vector(mode = "numeric", length = length(cell_types))
    names(p.temp) <- cell_types
    for (cell in cell_types) {
      dist <- dis[[trf]][[cell]]
      if (!mean(dist) == 0) {
        est <- fitdist((dist + 1), "lnorm", method = "mle")
        test <- ks.test(dist,
                        rlnorm(
                          n = length(dist),
                          meanlog = est$estimate[1],
                          sdlog = est$estimate[2]
                        ))
        p.temp[cell] <- test$p.value
      } else {
        p.temp[cell] <- 0
      }
    }
    pvals[[trf]] <- p.temp
    pos <- which(names(dis) == trf)
    if (pos %% 1000 == 1)
      message(pos)
  }
  saveRDS(pvals, "data/pvals_per_cell_ks_list.rds")
} else {
  pvals <- readRDS("data/pvals_per_cell_ks_list.rds")
}
pvals_df <- plyr::ldply(pvals)
rownames(pvals_df) <- pvals_df$.id
pvals_df$.id <- NULL
head(pvals_df)

pvals_df_bin <- data.frame(apply(pvals_df, 2, function(x) x>=.01)) #presence absence defined by ks test pvalue threshold
head(pvals_df_bin)

nrow(pvals_df_bin)

#>>intra vs. extracellular####
sum(apply(pvals_df_bin[, c("exosomes", "serum")], 1, any))
sum(apply(pvals_df_bin[, grep("cell", colnames(pvals_df_bin))], 1, any))
sum(pvals_df_bin$CD14.cells)
sum(pvals_df_bin$whole.blood)

#in DE
deseq <- read.csv("out/DESeq2_smrna_stroke_vs_con.csv")
pvals_df_bin_sig <- pvals_df_bin[rownames(pvals_df_bin) %in% deseq$name[deseq$padj < .05],]

sum(apply(pvals_df_bin_sig[, c("exosomes", "serum")], 1, any))
sum(apply(pvals_df_bin_sig[, grep("cell", colnames(pvals_df_bin))], 1, any))
sum(pvals_df_bin_sig$CD14.cells)
sum(pvals_df_bin_sig$whole.blood)

#>plot####
dat <- plotCounts(dds, "tRF-18-HR0VX6D2", intgroup = "cells", returnData = T)
dat$cells <- gsub(" ", ".", dat$cells, fixed = T)

#>>print example figure presence/absence####
#presence
cd14 <- dat[dat$cells == "CD14.cells",]
dist <- cd14$count
est <- fitdist((dist + 1), "lnorm", method = "mle")
test <- ks.test(dist,
                rlnorm(
                  n = length(dist),
                  meanlog = est$estimate[1],
                  sdlog = est$estimate[2]
                ))
ggplot(cd14) + geom_histogram(bins = 21, aes(x = count, y = ..density..), color = "grey") + 
  stat_function(fun = dlnorm, args = list(meanlog = est$estimate[1], sdlog = est$estimate[2]), 
                colour = "red") + 
  theme_light()
ggsave("img/presence_hr0v_cd14.svg")

#absence
ery <- dat[dat$cells == "CD235a.cells",]
dist <- ery$count
est <- fitdist((dist + 1), "lnorm", method = "mle")
test <- ks.test(dist,
                rlnorm(
                  n = length(dist),
                  meanlog = est$estimate[1],
                  sdlog = est$estimate[2]
                ))
pwr <- function(x) (.15*x^-.6)
ggplot(ery) + geom_histogram(bins = 22, 
  aes(x = count, y = ..density..), color = "grey") + 
  stat_function(fun = dlnorm, args = list(meanlog = est$estimate[1], sdlog = est$estimate[2]), 
                colour = "red") + 
  stat_function(fun = pwr, colour = "green") +
  theme_light()
ggsave("img/absence_hr0v_ery.svg")



# heatmap clustering####

pvals_df_num <- data.frame(apply(pvals_df_bin, 2, as.numeric))
rownames(pvals_df_num) <- rownames(pvals_df_bin)

cell_types
cell_types.name <- c("monocyte", "neutrophil", "B-cell", "erythrocyte", "T-helper",
                     "NK-cell", "T-toxic", "exosome", "serum", "whole.blood")

colnames(pvals_df_num) <- cell_types.name
mean(rowMeans(pvals_df_num)==0)
pvals_df_num <- pvals_df_num[rowMeans(pvals_df_num) != 0,]

dim(pvals_df_num)

# >ward.D2####
set.seed(1234)
method <- "ward.D2"
n <- 8
dev.off()
p <- pheatmap(
  as.matrix(pvals_df_num),
  color = colorRampPalette(brewer.pal(n = 3, name = "RdYlBu"))(2),
  clustering_method = method,
  legend = F,
  show_rownames = F,
  cutree_rows = n,
  # filename = paste0("img/heatmaps/", method, "-", n, ".pdf"),
  main = paste0(method, "-", n)
)

clust <- data.frame(cutree(as.hclust(p$tree_row), n))
colnames(clust) <- "cluster"
clust$cluster <- factor(clust$cluster)
#table out
clust0 <- clust
clust0$name <- rownames(clust0)
write.csv(clust0[, c("name", "cluster")], "out/tRF_cluster_table.csv", row.names = F)

pheatmap(as.matrix(pvals_df_num),
         color = colorRampPalette(brewer.pal(n = 3, name = "RdYlBu"))(2),
         annotation_colors = list(brewer.pal(n, "Dark2")),
         clustering_method = method,
         legend = T,
         show_rownames = F,
         annotation_row = clust,
         cutree_rows = n)

#clusters:
# 1 mono-NK
# 2 ubi
# 3 immune-minus-CD8
# 4 immune
# 5 blood-B
# 6 immune-minus-B
# 7 mono
# 8 B-T
levels(clust$cluster) <- c("mono-NK", "ubi", "immune-minus-CD8", "immune", "blood-B", "immune-minus-B", "mono", "B-T")

pheatmap(as.matrix(pvals_df_num),
         color = colorRampPalette(brewer.pal(n = 3, name = "RdYlBu"))(2),
         annotation_colors = list(brewer.pal(n, "Dark2")),
         clustering_method = method,
         legend = T,
         show_rownames = F,
         annotation_row = clust,
         cutree_rows = n,
         filename = "img/trf_heatmap_cell_type.pdf"); dev.off()


#t-SNE of PRESENCE/ABSENCE####
deseq <- read.csv("out/DESeq2_smrna_stroke_vs_con.csv")

#>cluster from pheatmap####
tsne_trf <- Rtsne(as.matrix(pvals_df_num), check_duplicates = F, perplexity = 50)
tsne_trf_plot <- data.frame(tsne_trf$Y)
colnames(tsne_trf_plot) <- c("x", "y")
rownames(tsne_trf_plot) <- rownames(pvals_df_num)
tsne_trf_plot$shape <- factor(rownames(tsne_trf_plot) %in% deseq$name)
tsne_trf_plot$countChange <- deseq$countChange[match(rownames(tsne_trf_plot), deseq$name)]

#>>ward.D2_8####
col <- clust$cluster[match(rownames(pvals_df_num), rownames(clust))]
tsne_trf_plot$ward.D2_8 <- factor(col)

#complete
ggplot(tsne_trf_plot, aes(x, y, col = ward.D2_8)) + geom_point(size = 1) +
  scale_color_brewer(palette = "Dark2") + theme_minimal()
ggsave("img/tsne_trf.svg", width = 8, height = 5)

#de
ggplot(tsne_trf_plot, aes(x, y)) + geom_point(size = 1, col = "grey80") +
  theme_minimal() +
  geom_point(aes(x, y, col = ward.D2_8, size = countChange), data =
               tsne_trf_plot[tsne_trf_plot$shape == T,]) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(breaks = c(10, 30, 100, 300, 1000, 3000))
ggsave("img/tsne_trf_DE.svg", width = 8, height = 5)

#cholinergic
trf_counts <- read.csv("out/cholinotrf_counts.csv")
tsne_trf_plot$name <- rownames(tsne_trf_plot)
ggplot(tsne_trf_plot, aes(x, y)) + geom_point(size = 1, col = "grey80") +
  theme_minimal() +
  geom_point(aes(x, y, col = ward.D2_8, size = countChange), data =
               tsne_trf_plot[rownames(tsne_trf_plot) %in% trf_counts$MINTplate[trf_counts$n>4],]) +
  geom_text_repel(aes(x, y, label = name, size = countChange), data = tsne_trf_plot[rownames(tsne_trf_plot) %in% trf_counts$MINTplate[trf_counts$n>4],]) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(breaks = c(10, 30, 100, 300, 1000, 3000))
ggsave("img/tsne_trf_cholino.svg", width = 8, height = 5)
