setwd("~/GitHub/stroke-trf")
rm(list = ls())

pacman::p_load(Rtsne,
               reshape2,
               RColorBrewer,
               ggplot2)

options(stringsAsFactors = F)

cholinergic_genes <- readRDS("data/cholinergic_genes.rds")

#tissues from marbach et al that are equivalent to juzenas et al tissues
tissues <- c("CD4POS_T_CELLS", "CD8POS_T_CELLS", "CD14POS_MONOCYTES", "NATURAL_KILLER_CELLS", "NEUTROPHILS", "CD19POS_B_CELLS", "WHOLE_BLOOD")

# table with cholinergic expression in each tissue####
# tissue files have been transformed from raw marbach et al 2016 data into .rds format
# expression is calculated by summation of transcription factor activity towards each gene per tissue
files <- paste0("data/", list.files("data/"))
ch_exp <- cholinergic_genes
allgenes <- character()
for(tis in tissues){
  idx <- which(tissues == tis)
  message(idx)
  file <- files[grep(tis, files)]
  tar <- readRDS(file)
  tar_sum <- aggregate(tar$r.tfa, by = list(g.name = tar$g.name, g.ensg = tar$g.ensg), "sum")
  tar_sum <- tidyr::separate_rows(tar_sum, g.ensg, sep = ", ")
  
  #cholinergic
  col <- tar_sum$x[match(ch_exp$ensg, tar_sum$g.ensg)]
  col[is.na(col)] <- 0
  ch_exp <- cbind(ch_exp, col)
  colnames(ch_exp)[idx+4] <- tis
  
  #prepare complete
  allgenes <- unique(c(allgenes, tar_sum$g.ensg))
}
head(ch_exp)

all_exp <- data.frame(ensg = allgenes)
all_exp$name <- NA
for(tis in tissues){
  idx <- which(tissues == tis)
  message(idx)
  file <- files[grep(tis, files)]
  tar <- readRDS(file)
  tar_sum <- aggregate(tar$r.tfa, by = list(g.name = tar$g.name, g.ensg = tar$g.ensg), "sum")
  tar_sum <- tidyr::separate_rows(tar_sum, g.ensg, sep = ", ")
  
  #cholinergic
  col <- tar_sum$x[match(all_exp$ensg, tar_sum$g.ensg)]
  col[is.na(col)] <- 0
  all_exp <- cbind(all_exp, col)
  namidx <- which(all_exp$ensg %in% tar_sum$g.ensg)
  all_exp$name[namidx] <- tar_sum$g.name[match(all_exp$ensg[namidx], tar_sum$g.ensg)]
  colnames(all_exp)[idx+2] <- tis
  
  #prepare complete
  allgenes <- unique(c(allgenes, tar_sum$g.ensg))
}
head(all_exp)
nrow(all_exp)

#t-SNE of expression####
which(rowMeans(all_exp[3:ncol(all_exp)]) == 0)
#>tissues akin to juzenas? merge? select?####
tissues
datExp <- all_exp[, 3:ncol(all_exp)]
dim(datExp)
rownames(datExp) <- all_exp$ensg

#color by highest expression (derived from transcriptional activity towards each gene)####
#including whole blood
datExp0 <- datExp#[1:2000,]
head(datExp0)

tsne_tf <-
  Rtsne(as.matrix(datExp0),
        check_duplicates = F,
        perplexity = 43)
tsne_tf_plot <- data.frame(tsne_tf$Y)
colnames(tsne_tf_plot) <- c("x", "y")
rownames(tsne_tf_plot) <- rownames(datExp0)
tsne_tf_plot$name <- rownames(datExp0)
tsne_tf_plot$name[tsne_tf_plot$name %in% cholinergic_genes$ensg] <-
  cholinergic_genes$gene_symbol[match(tsne_tf_plot$name[tsne_tf_plot$name %in% cholinergic_genes$ensg], cholinergic_genes$ensg)]

# find shared expression
find.max <- function(row) {
  return(names(row)[which.max(row)])
}

tsne_tf_plot$cluster <- apply(datExp0, 1, find.max)

ggplot(tsne_tf_plot, aes(x, y, col = cluster)) + geom_point(size = 2) +
  theme_minimal()
ggsave(paste0("img/tsne_tfs_all.pdf"))

ggplot(tsne_tf_plot, aes(x, y)) + geom_point(size = 2, col = "grey80") +
  geom_point(aes(x, y, col = cluster), data = tsne_tf_plot[tsne_tf_plot$name %in% cholinergic_genes$gene_symbol[cholinergic_genes$group %in% c("core", "receptor")],]) +
  geom_text(aes(x, y, label = name), data = tsne_tf_plot[tsne_tf_plot$name %in% cholinergic_genes$gene_symbol[cholinergic_genes$group %in% c("core", "receptor")],]) +
  theme_minimal()
ggsave(paste0("img/tsne_tfs_cholinergic.pdf"))


#without whole blood####
datExp0 <- datExp[,colnames(datExp) != "WHOLE_BLOOD"]
datExp0 <- datExp0[which(rowMeans(datExp0) > 0), ]
dim(datExp0)

set.seed(1)
tsne_tf <-
  Rtsne(as.matrix(datExp0),
        check_duplicates = F,
        perplexity = 44)
tsne_tf_plot <- data.frame(tsne_tf$Y)
colnames(tsne_tf_plot) <- c("x", "y")
rownames(tsne_tf_plot) <- rownames(datExp0)
tsne_tf_plot$name <- rownames(datExp0)
tsne_tf_plot$name[tsne_tf_plot$name %in% cholinergic_genes$ensg] <-
  cholinergic_genes$gene_symbol[match(tsne_tf_plot$name[tsne_tf_plot$name %in% cholinergic_genes$ensg], cholinergic_genes$ensg)]

# find shared expression
find.max <- function(row) {
  return(names(row)[which.max(row)])
}

tsne_tf_plot$cluster <- apply(datExp0, 1, find.max)

ggplot(tsne_tf_plot, aes(x, y, col = cluster)) + geom_point(size = 1) +
  theme_minimal()
ggsave(
  paste0("img/tsne_tfs_all_no_whole.svg"),
  width = 12,
  height = 10
)

ggplot(tsne_tf_plot, aes(x, y)) + geom_point(size = 1, col = "grey80") +
  geom_point(aes(x, y, col = cluster, size = 3), data = tsne_tf_plot[tsne_tf_plot$name %in% cholinergic_genes$gene_symbol[cholinergic_genes$group %in% c("core", "receptor")],]) +
  geom_text(aes(x, y, label = name), data = tsne_tf_plot[tsne_tf_plot$name %in% cholinergic_genes$gene_symbol[cholinergic_genes$group %in% c("core", "receptor")],]) +
  theme_minimal()
ggsave(paste0("img/tsne_tfs_cholinergic_no_whole.svg"), width = 12, height = 10)
