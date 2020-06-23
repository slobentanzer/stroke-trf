setwd("~/GitHub/stroke-trf")
rm(list = ls())

library(pheatmap)
library(fitdistrplus)
library(DESeq2)
library(ggplot2)
library(Rtsne)
library(RColorBrewer)
library(ggrepel)

options(stringsAsFactors = FALSE);

#Read in
meta <- read.csv("data/raw/juzenas_sample_meta.csv", sep = ",", row.names = 1)
meta$cells <- factor(meta$cells)
head(meta)
rownames(meta) <- meta$file
meta$batch <- factor(meta$batch)

counts <- read.csv("data/raw/miRNA_reads_miRExpress_juzenas.csv", sep = ",", row.names = 1)
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
if(!file.exists("data/counts_per_cell_list_mir.rds")){
  dis <- vector(mode = "list", length = nrow(cts))
  names(dis) <- rownames(cts)
  for (mir in rownames(counts)) {
    row <- counts[mir,]
    expr <- logical()
    lis <- vector(mode = "list", length = length(unique(meta$cells)))
    names(lis) <- unique(meta$cells)
    for (cell in unique(meta$cells)) {
      lis[[cell]] <- as.numeric(row[meta$file[meta$cells == cell]])
    }
    dis[[mir]] <- lis
  }
  saveRDS(dis, "data/counts_per_cell_list_mir.rds")
} else {
  dis <- readRDS("data/counts_per_cell_list_mir.rds")
}

#>implement statistical test####
cell_types <- levels(meta$cells)

#kolmogorov-smirnov over log-normal distribution of counts with estimated mean and sd
if(!file.exists("data/pvals_per_cell_ks_list_mir.rds")){
  pvals <- vector(mode = "list", length = length(dis))
  names(pvals) <- names(dis)
  
  
  for(mir in names(dis)) {
    p.temp <- vector(mode = "numeric", length = length(cell_types))
    names(p.temp) <- cell_types
    for (cell in cell_types) {
      dist <- dis[[mir]][[cell]]
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
    pvals[[mir]] <- p.temp
    pos <- which(names(dis) == mir)
    if (pos %% 1000 == 1)
      message(pos)
  }
  saveRDS(pvals, "data/pvals_per_cell_ks_list_mir.rds")
} else {
  pvals <- readRDS("data/pvals_per_cell_ks_list_mir.rds")
}
pvals_df <- plyr::ldply(pvals)
rownames(pvals_df) <- pvals_df$.id
pvals_df$.id <- NULL
head(pvals_df)

pvals_df_bin <- data.frame(apply(pvals_df, 2, function(x) x>=.01)) #presence absence defined by ks test pvalue threshold
head(pvals_df_bin)

pvals_df_num <- data.frame(apply(pvals_df_bin, 2, as.numeric))
rownames(pvals_df_num) <- rownames(pvals_df_bin)

cell_types
cell_types.name <- c("monocyte", "neutrophil", "B-cell", "erythrocyte", "T-helper", 
                     "NK-cell", "T-toxic", "exosome", "serum", "whole.blood")

colnames(pvals_df_num) <- cell_types.name
mean(rowMeans(pvals_df_num)==0)
pvals_df_num <- pvals_df_num[rowMeans(pvals_df_num) != 0,]

dim(pvals_df_num)

# heatmap clustering####

# ward.D 8
set.seed(1234)
method <- "ward.D"
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

pheatmap(as.matrix(pvals_df_num),
         color = colorRampPalette(brewer.pal(n = 3, name = "RdYlBu"))(2),
         annotation_colors = list(brewer.pal(n, "Dark2")),
         clustering_method = method,
         legend = T,
         show_rownames = F,
         annotation_row = clust,
         cutree_rows = n)

#clusters:
# 1 immune-neutro-whole
# 2 ubi
# 3 immune-neutro-whole-checker
# 4 ubi-checker
# 5 immune
# 6 mono
# 7 ubi-checker-less
# 8 whole-blood
levels(clust$cluster) <- c("immune-neutro-whole", "ubi", "immune-neutro-whole-checker", "ubi-checker", "immune", "mono", "ubi-checker-less", "whole-blood")

pheatmap(as.matrix(pvals_df_num),
         color = colorRampPalette(brewer.pal(n = 3, name = "RdYlBu"))(2),
         annotation_colors = list(brewer.pal(n, "Dark2")),
         clustering_method = method,
         legend = T,
         show_rownames = F,
         annotation_row = clust,
         cutree_rows = n,
         filename = "img/mir_heatmap_cell_type.pdf"); dev.off()

#t-SNE of PRESENCE/ABSENCE####
deseq <- read.csv("out/DESeq2_smrna_stroke_vs_con.csv")

#>cluster from pheatmap####
tsne_mir <- Rtsne(as.matrix(pvals_df_num), check_duplicates = F, perplexity = 50)
tsne_mir_plot <- data.frame(tsne_mir$Y)
colnames(tsne_mir_plot) <- c("x", "y")
rownames(tsne_mir_plot) <- rownames(pvals_df_num)
tsne_mir_plot$shape <- factor(rownames(tsne_mir_plot) %in% deseq$name)
tsne_mir_plot$countChange <- abs(deseq$countChange[match(rownames(tsne_mir_plot), deseq$name)])
tsne_mir_plot$log_cc <- log10(tsne_mir_plot$countChange+2)

#>>ward.D_8####
col <- clust$cluster[match(rownames(pvals_df_num), rownames(clust))]
tsne_mir_plot$ward.D_8 <- factor(col)

#complete
ggplot(tsne_mir_plot, aes(x, y, col = ward.D_8)) + geom_point(size = 1) +
  scale_color_brewer(palette = "Dark2") + theme_minimal()
# geom_point(aes(x, y), col = "white", pch = 4, data =
#              tsne_mir_plot[tsne_mir_plot$shape == T & !rownames(tsne_mir_plot) %in% top20,]) +
# geom_point(aes(x, y), col = "black", pch = 4, data = tsne_mir_plot[rownames(tsne_mir_plot) %in% top20,])
ggsave("img/tsne_mir.svg", width = 8, height = 5)

#de
ggplot(tsne_mir_plot, aes(x, y)) + geom_point(size = 1, col = "grey80") +
  theme_minimal() +
  geom_point(aes(x, y, col = ward.D_8, size = log_cc), data =
               tsne_mir_plot[tsne_mir_plot$shape == T,]) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(breaks = 1:8, labels = 10^(1:8), name = "count change", range = c(.4, 3))
ggsave("img/tsne_mir_DE.svg", width = 8, height = 5)

#cholinergic
mir_counts <- read.csv("out/cholinomiR_counts.csv")
tsne_mir_plot$name <- rownames(tsne_mir_plot)
ggplot(tsne_mir_plot, aes(x, y)) + geom_point(size = 1, col = "grey80") +
  theme_minimal() +
  geom_point(aes(x, y, col = ward.D_8, size = log_cc), data =
               tsne_mir_plot[rownames(tsne_mir_plot) %in% mir_counts$m.name[mir_counts$n>4],]) +
  geom_text_repel(aes(x, y, label = name, size = log_cc), data = tsne_mir_plot[rownames(tsne_mir_plot) %in% mir_counts$m.name[mir_counts$n>4],]) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_continuous(breaks = 1:4, labels = 10^(1:4), name = "count change", range = c(.4, 3))
ggsave("img/tsne_mir_cholino.svg", width = 8, height = 5)
