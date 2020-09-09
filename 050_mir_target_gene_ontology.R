# NETWORK GENERATION WAS PERFORMED USING MIRNEO (INHOUSE DATABASE, see Lobentanzer et al. 2019, Cell Reports) #
# for execution of Neo4j queries, a local version of the database must be run #
# temporary save files for the network queries are supplied for reproduction without access to the database #

setwd("~/GitHub/stroke-trf")
rm(list = ls())

db_running <- F #access to local miRNeo database?

pacman::p_load(RNeo4j, BiocParallel, gsoap, RColorBrewer, ggplot2, tidyverse, janitor, car, lsmeans)

if(db_running){
  graph <- startGraph("http://localhost:7474/db/data/")
  query <- "MATCH (m:MIR {species: 'HSA'}) RETURN m.name"
  allmirs <- cypher(graph, query)$m.name
  query <- "MATCH (g:GENE {species: 'HSA'}) RETURN g.ensg"
  allgenes <- cypher(graph, query)$g.ensg
} else {
  load("data/complete_mir_gene_list.RData")
}
ensg2symbol <- readRDS(file = "data/ensg2symbol_all.rds")

source("GO.R")

# smRNA GO ANALYSIS ####
de_smrna <- read.csv(file = "out/DESeq2_smrna_stroke_vs_con.csv")
de_mirs <- de_smrna[-grep("tRF", de_smrna$name),]
# write.csv(de_mirs, file = "out/DESeq2_mirna_stroke_vs_con.csv")

de_genes <- read.csv(file = "out/DESeq2_gene_stroke_vs_con.csv")

# >functions ####
get.targeting <- function(mirs, genes = NULL) {
  if(is.null(genes)){
    targets <- cypher(graph, "MATCH (m:MIR)-[r:ALGOSUM]-(g:GENE {species: 'HSA'}) 
                  WHERE m.name in {mirnames}
                  AND r.rating > 5
                  RETURN m.name, r.rating, g.name, g.id, g.ensg",
                      mirnames = mirs)
  } else {
    targets <- cypher(graph, "MATCH (m:MIR)-[r:ALGOSUM]-(g:GENE {species: 'HSA'}) 
                  WHERE m.name in {mirnames} AND g.name in {genenames}
                  AND r.rating > 5
                  RETURN m.name, r.rating, g.name, g.id, g.ensg",
                      mirnames = mirs, genenames = genes)
  }
}

permute.targeting <- function(num, mirlist, genelist = NULL, path) {
  mirs_perm <- sample(allmirs, length(mirlist))
  if(is.null(genelist)){
    targets_perm <- get.targeting(mirs_perm)
  } else {
    targets_perm <- get.targeting(mirs_perm, genes = genelist)
  }
  tar_count_perm <- dplyr::count(targets_perm, g.ensg, sort = T)
  tar_count_perm <- tar_count_perm[!is.na(tar_count_perm$g.ensg),]
  if(!dir.exists(path))
    dir.create(path, recursive = T)
  write.table(tar_count_perm, row.names = F, sep = ";", quote = F, file = paste0(path, num, ".txt"))
}

get.pval <- function(row){
  row <- as.numeric(row)
  tr.val <- row[1]
  dist <- sort(row[-1])
  dist <- cbind(null = dist, data.frame(prob=1-abs(2*seq(1:length(dist))/length(dist)-1)))
  pval <- dist$prob[findInterval(tr.val, dist$null)]
  if(length(pval) == 1){
    return(pval)
  } else {
    return(0)
  }
}


#gsoap - get single significant genes
get.genes.in.go.id <- function(id, GOdata, test = NULL, module = NULL) {
  if(!is.null(test)) {
    if (!is.null(module)) {
      sigGenes <- test[test$ensg %in% unlist(genesInTerm(GOdata, id)) &
                         test$ffl_mod == module & test$padj < .05,]
    } else {
      sigGenes <- test[test$ensg %in% unlist(genesInTerm(GOdata, id)) &
                         test$padj < .05,]
    }
    return(paste(sigGenes$symbol, collapse = "/"))
  } else {
    sigGenes <- unlist(genesInTerm(GOdata, id))
    return(paste(sigGenes, collapse = "/"))
  }
}

ShapiroTest <- function(x){
  z <- length(x)
  if(z > 2 && z < 5001){
    y <- shapiro.test(x)
    return (c(w=unname(y$statistic), p=y$p.value, n=z))
  }
  return (c(w=NA, p=NA, n=z))
}

OneWayAnova <- function(data, dv, iv, conf = 0.95){
  
  alpha = 1 - conf
  ndigits = 10
  
  df <- data
  
  # empty result object
  r <- {}
  dfFormula <- as.formula(paste(dv, iv, sep=" ~ "))
  df.res <- lm(dfFormula, data = df)
  r$aov <- do.call(rbind.data.frame, summary(aov(df.res)))
  
  r$aov <- cbind(Source=rownames(r$aov), r$aov)
  
  suppressWarnings(suppressMessages(library(lsmeans)))
  
  t<- suppressMessages(pairs(lsmeans(df.res, iv)))
  
  ts <- summary(t)
  tc <- confint(t, level = conf)
  r$t <- merge(ts, tc)  
  
  r$s <- aggregate(formula=dfFormula,
                   data=df,
                   FUN=ShapiroTest)

  return(r)
}

# >targeting ####
nperm <- 1e4

if(db_running){
  # >>> positive reg###
  mirs_pos <- de_mirs$name[de_mirs$log2FoldChange > 0 & de_mirs$padj < .05]
  length(mirs_pos)
  targets_pos <- get.targeting(mirs_pos)
  # >>> negative reg###
  mirs_neg <- de_mirs$name[de_mirs$log2FoldChange < 0 & de_mirs$padj < .05]
  length(mirs_neg)
  targets_neg <- get.targeting(mirs_neg)
  # >>> negative with threshold###
  mirs_neg_logfc2 <- de_mirs$name[de_mirs$log2FoldChange < -2 & de_mirs$padj < .05]
  length(mirs_neg_logfc2)
  targets_neg_logfc2 <- get.targeting(mirs_neg_logfc2)
  save(targets_pos, targets_neg, targets_neg_logfc2, file = "data/mir_targets_graph.RData")
} else {
  load("data/mir_targets_graph.RData")
}
  
# >>> positive reg####
head(targets_pos)
nrow(targets_pos)
length(unique(targets_pos$g.name)) #number of unique targets_pos
tar_count_pos <- dplyr::count(targets_pos, g.ensg, sort = T) #sort by number of mirs targeting each gene
tar_count_pos <- tar_count_pos[!is.na(tar_count_pos$g.ensg),]
head(tar_count_pos, 20)
range(tar_count_pos$n)
hist(tar_count_pos$n)
mean(tar_count_pos$n)
sd(tar_count_pos$n)

# targets per mir 
mir_count_pos <- dplyr::count(targets_pos, m.name, sort = T)
head(mir_count_pos, 20)

# >>> negative reg####
head(targets_neg)
nrow(targets_neg)
length(unique(targets_neg$g.name)) #number of unique targets_neg
tar_count_neg <- dplyr::count(targets_neg, g.ensg, sort = T) #sort by number of mirs targeting each gene
tar_count_neg <- tar_count_neg[!is.na(tar_count_neg$g.ensg),]
head(tar_count_neg, 20)
range(tar_count_neg$n)
hist(tar_count_neg$n)
mean(tar_count_neg$n)
sd(tar_count_neg$n)

# targets per mir 
mir_count_neg <- dplyr::count(targets_neg, m.name, sort = T)
head(mir_count_neg, 20)

# >>>> diverging number of targets in pos and neg regulated miRs ####
hist(mir_count_neg$n)
hist(mir_count_pos$n)
hist(log10(mir_count_neg$n))
hist(log10(mir_count_pos$n))
mean(mir_count_pos$n)
mean(mir_count_neg$n)
median(mir_count_pos$n)
median(mir_count_neg$n)

mir_count_neg$group <- "Negative"
mir_count_pos$group <- "Positive"

mir_count <- rbind(mir_count_neg, mir_count_pos)
mir_count$nlog <- log10(mir_count$n)
OneWayAnova(mir_count, dv = "n", iv = "group")
OneWayAnova(mir_count, dv = "nlog", iv = "group")

ggplot(mir_count, aes(group, log10(n), fill = group)) + geom_boxplot() +
  scale_y_continuous(breaks = 0:4, labels = c(1, 10, 100, 1000, 10000), name = "Number of Targets") +
  scale_x_discrete(name = "Direction of Regulation") + theme(legend.position = "none") +
  theme_minimal()
ggsave("img/mir_target_numbers_by_direction_of_regulation.pdf", width = 6, height = 4)

# >>>> positive: permute inside graph ####
if(db_running){
  # to remove genes that are not enriched in targeting
  dir.create("data/permutation/mirs_pos/", recursive = T)
  done <-
    as.numeric(gsub(".txt", "", fixed = T, list.files("data/permutation/mirs_pos/")))
  iter <- 1:nperm
  remain <- iter[!iter %in% done]
  
  bplapply(remain,
           permute.targeting,
           mirlist = mirs_pos,
           path = "data/permutation/mirs_pos/", BPPARAM = MulticoreParam(workers = 7))
  
  # >>>>> positive: read permutation ####
  gene_table_pos <-
    data.frame(row.names = allgenes[!is.na(allgenes) &
                                      !duplicated(allgenes)])
  gene_table_pos$res <-
    tar_count_pos$n[match(rownames(gene_table_pos), tar_count_pos$g.ensg)]
  gene_table_pos <- subset(gene_table_pos,!is.na(gene_table_pos$res))
  
  files <- list.files("data/permutation/mirs_pos", full.names = T)
  for (file in files) {
    df <- data.table::fread(file, sep = ";")
    vec <- df$n[match(rownames(gene_table_pos), df$g.ensg)]
    gene_table_pos[, (ncol(gene_table_pos) + 1)] <- vec
  }
  gene_table_pos[is.na(gene_table_pos)] <- 0
  gene_table_pos <-
    gene_table_pos[order(gene_table_pos$res, decreasing = T), ]
  
  #row-wise p-values
  gene_table_pos$pval <- unlist(apply(gene_table_pos, 1, get.pval))
  gene_table_pos$mean <-
    rowMeans(gene_table_pos[,!colnames(gene_table_pos) %in% c("res", "pval")])
  sig_table_pos <-
    gene_table_pos[gene_table_pos$pval < .05 &
                     gene_table_pos$res > gene_table_pos$mean, c("res", "pval", "mean")]
  nrow(sig_table_pos)
  
  # fill with nonsig genes for background
  if (nrow(sig_table_pos) < 5000) {
    set.seed(123)
    nams <-
      sample(rownames(gene_table_pos)[!rownames(gene_table_pos) %in% rownames(sig_table_pos)], 5000 -
               nrow(sig_table_pos))
    sig_table_pos <-
      rbind(sig_table_pos, gene_table_pos[nams, c("res", "pval", "mean")])
    sig_table_pos$res[rownames(sig_table_pos) %in% nams] <- 0
  }
  
  go_pos_top200 <-
    get.GO(sig_table_pos$res,
           rownames(sig_table_pos),
           quant = .96,
           top = 5000)
  go_pos_top200$result
  
  go_pos_by_pval <-
    get.GO(sig_table_pos$pval,
           rownames(sig_table_pos),
           pval = T,
           top = 5000)
  go_pos_by_pval$result
  
  # >>>> negative: permute inside graph ####
  # to remove genes that are not enriched in targeting
  dir.create("data/permutation/mirs_neg/", recursive = T)
  done <-
    as.numeric(gsub(".txt", "", fixed = T, list.files("data/permutation/mirs_neg/")))
  iter <- 1:nperm
  remain <- iter[!iter %in% done]
  
  bplapply(remain,
           permute.targeting,
           mirlist = mirs_neg,
           path = "data/permutation/mirs_neg/", BPPARAM = MulticoreParam(workers = 7))
  
  # >>>>> negative: read permutation ####
  gene_table_neg <-
    data.frame(row.names = allgenes[!is.na(allgenes) &
                                      !duplicated(allgenes)])
  gene_table_neg$res <-
    tar_count_neg$n[match(rownames(gene_table_neg), tar_count_neg$g.ensg)]
  gene_table_neg <- subset(gene_table_neg,!is.na(gene_table_neg$res))
  
  files <- list.files("data/permutation/mirs_neg", full.names = T)
  for (file in files) {
    df <- data.table::fread(file, sep = ";")
    vec <- df$n[match(rownames(gene_table_neg), df$g.ensg)]
    gene_table_neg[, (ncol(gene_table_neg) + 1)] <- vec
  }
  gene_table_neg[is.na(gene_table_neg)] <- 0
  gene_table_neg <-
    gene_table_neg[order(gene_table_neg$res, decreasing = T), ]
  
  #row-wise p-values
  gene_table_neg$pval <- unlist(apply(gene_table_neg, 1, get.pval))
  gene_table_neg$mean <-
    rowMeans(gene_table_neg[,!colnames(gene_table_neg) %in% c("res", "pval")])
  sig_table_neg <-
    gene_table_neg[gene_table_neg$pval < .05 &
                     gene_table_neg$res > gene_table_neg$mean, c("res", "pval", "mean")]
  nrow(sig_table_neg)
  
  # fill with nonsig genes for background
  if (nrow(sig_table_neg) < 5000) {
    set.seed(123)
    nams <-
      sample(rownames(gene_table_neg)[!rownames(gene_table_neg) %in% rownames(sig_table_neg)], 5000 -
               nrow(sig_table_neg))
    sig_table_neg <-
      rbind(sig_table_neg, gene_table_neg[nams, c("res", "pval", "mean")])
    sig_table_neg$res[rownames(sig_table_neg) %in% nams] <- 0
  }
  
  go_neg_top200 <-
    get.GO(sig_table_neg$res,
           rownames(sig_table_neg),
           quant = .96,
           top = 5000)
  go_neg_top200$result
  
  # more than 5000 genes with pval < .05 -> threshold?
  quantile(sig_table_neg$pval, seq(0, 1, .1))
  table(sig_table_neg$pval)
  
  # >> logFC threshold for negatively reg. mirs ####
  # >>> true value: negative reg####
  head(targets_neg_logfc2)
  nrow(targets_neg_logfc2)
  length(unique(targets_neg_logfc2$g.name)) #number of unique targets
  tar_count_neg_logfc2 <- dplyr::count(targets_neg_logfc2, g.ensg, sort = T) #sort by number of mirs targeting each gene
  tar_count_neg_logfc2 <- tar_count_neg_logfc2[!is.na(tar_count_neg_logfc2$g.ensg),]
  head(tar_count_neg_logfc2, 20)
  range(tar_count_neg_logfc2$n)
  hist(tar_count_neg_logfc2$n)
  mean(tar_count_neg_logfc2$n)
  sd(tar_count_neg_logfc2$n)
  
  # >>> negative: permute inside graph ####
  # to remove genes that are not enriched in targeting
  dir.create("data/permutation/mirs_neg_logfc2/", recursive = T)
  done <- as.numeric(gsub(".txt", "", fixed = T, list.files("data/permutation/mirs_neg_logfc2/")))
  iter <- 1:nperm
  remain <- iter[!iter %in% done]
  
  bplapply(remain, permute.targeting, mirlist = mirs_neg_logfc2, path = "data/permutation/mirs_neg_logfc2/", BPPARAM = MulticoreParam(workers = 7))
  
  # >>>>> negative: read permutation ####
  gene_table_neg_logfc2 <- data.frame(row.names = allgenes[!is.na(allgenes) & !duplicated(allgenes)])
  gene_table_neg_logfc2$res <- tar_count_neg_logfc2$n[match(rownames(gene_table_neg_logfc2), tar_count_neg_logfc2$g.ensg)]
  gene_table_neg_logfc2 <- subset(gene_table_neg_logfc2, !is.na(gene_table_neg_logfc2$res))
  
  files <- list.files("data/permutation/mirs_neg_logfc2", full.names = T)
  for(file in files){
    df <- data.table::fread(file, sep = ";")
    vec <- df$n[match(rownames(gene_table_neg_logfc2), df$g.ensg)]
    gene_table_neg_logfc2[, (ncol(gene_table_neg_logfc2)+1)] <- vec
  }
  gene_table_neg_logfc2[is.na(gene_table_neg_logfc2)] <- 0
  gene_table_neg_logfc2 <- gene_table_neg_logfc2[order(gene_table_neg_logfc2$res, decreasing = T),]
  
  #row-wise p-values
  gene_table_neg_logfc2$pval <- unlist(apply(gene_table_neg_logfc2, 1, get.pval))
  gene_table_neg_logfc2$mean <- rowMeans(gene_table_neg_logfc2[, !colnames(gene_table_neg_logfc2) %in% c("res", "pval")])
  sig_table_neg_logfc2 <- gene_table_neg_logfc2[gene_table_neg_logfc2$pval < .05 & gene_table_neg_logfc2$res > gene_table_neg_logfc2$mean, c("res", "pval", "mean")]
  nrow(sig_table_neg_logfc2)
  
  # fill with nonsig genes for background
  if(nrow(sig_table_neg_logfc2)<5000){
    set.seed(123)
    nams <- sample(rownames(gene_table_neg_logfc2)[!rownames(gene_table_neg_logfc2) %in% rownames(sig_table_neg_logfc2)], 5000-nrow(sig_table_neg_logfc2))
    sig_table_neg_logfc2 <- rbind(sig_table_neg_logfc2, gene_table_neg_logfc2[nams, c("res", "pval", "mean")])
    sig_table_neg_logfc2$res[rownames(sig_table_neg_logfc2) %in% nams] <- 0
  }
  
  go_neg_logfc2_by_pval <- get.GO(sig_table_neg_logfc2$pval, rownames(sig_table_neg_logfc2), pval = T, p.threshold = .0001, top = 5000)
  go_neg_logfc2_by_pval$result
  
  #save temp####
  save(go_pos_top200, go_pos_by_pval,
       go_neg_top200, go_neg_logfc2_by_pval,
       file = "data/mir_target_godata.RData")
} else {
  #load
  load("data/mir_target_godata.RData")
}

#GSOAP VISUALISATION####
#compile
df1 <- go_pos_top200$result
df1$genes <- unlist(lapply(df1$GO.ID, get.genes.in.go.id, GOdata = go_pos_top200$GOdata))
df1$weight <- as.numeric(df1$weight)
rownames(df1) <- df1$Term
head(df1)

df2 <- go_pos_by_pval$result
df2$genes <- unlist(lapply(df2$GO.ID, get.genes.in.go.id, GOdata = go_pos_by_pval$GOdata))
df2$weight <- as.numeric(df2$weight)
rownames(df2) <- df2$Term
head(df2)

df1$group <- "top200"
rownames(df1) <- paste(rownames(df1), df1$group, sep = "_")
df2$group <- "pval"
rownames(df2) <- paste(rownames(df2), df2$group, sep = "_")

df_pos <- rbind(df1, df2)
head(df_pos)

df3 <- go_neg_top200$result
df3$genes <- unlist(lapply(df3$GO.ID, get.genes.in.go.id, GOdata = go_neg_top200$GOdata))
df3$weight <- as.numeric(df3$weight)
rownames(df3) <- df3$Term
head(df3)

df4 <- go_neg_logfc2_by_pval$result
df4$genes <- unlist(lapply(df4$GO.ID, get.genes.in.go.id, GOdata = go_neg_logfc2_by_pval$GOdata))
df4$weight <- as.numeric(df4$weight)
idx <- which(duplicated(df4$Term) | duplicated(df4$Term, fromLast = T))
df4$Term[idx] <- paste0(df4$Term[idx], df4$GO.ID[idx])
rownames(df4) <- df4$Term
head(df4)

df3$group <- "top200"
rownames(df3) <- paste(rownames(df3), df3$group, sep = "_")
df4$group <- "pval"
rownames(df4) <- paste(rownames(df4), df4$group, sep = "_")

df_neg <- rbind(df3, df4)
head(df_neg)

# >compare between categories ####
df_pos$group <- paste(df_pos$group, "pos", sep = "_")
df_neg$group <- paste(df_neg$group, "neg", sep = "_")

df <- rbind(df_pos, df_neg)

rownames(df) <- paste0(df$GO.ID, ": ", df$Term, "_", df$group)
head(df)

#>gsoap layout####
set.seed(123)
l <- gsoap_layout(df,
                  "genes",
                  "weight",
                  projection = "tsne",
                  tsne.perplexity = 15)
quantile(l$significance)
idx <- seq(nrow(l))
l$group <-
  rownames(l) %>% strsplit("_") %>% lapply("[", c(2, 3)) %>% lapply(paste, collapse = "_") %>% unlist()
l$group <-
  factor(l$group,
         levels = c("pval_pos", "top200_pos", "pval_neg", "top200_neg"))
rownames(l) <-
  rownames(l) %>% gsub("_pval_pos", " ", .) %>% gsub("_pval_neg", ",", .) %>%
  gsub("_top200_pos", "", .) %>% gsub("_top200_neg", ".", ., fixed = T)
gsoap_plot(
  l,
  as.alpha = 'significance',
  as.color = 'group',
  which.labels = idx,
  # size.guide.loc = c(1., 1.),
  label.fontsize = 8,
  repel.direction = "both"
) + scale_color_manual(values = brewer.pal(6, "Paired")[3:6]) + scale_fill_manual(values = brewer.pal(6, "Paired")[3:6])
ggsave(
  "img/gsoap_gene_targets_pos_vs_neg_reg_mirnas.pdf",
  width = 30,
  height = 20
)
#without labels
gsoap_plot(
  l,
  as.alpha = 'significance',
  as.color = 'group',
  # which.labels = idx,
  # size.guide.loc = c(1., 1.),
  # label.fontsize = 8,
  # repel.direction = "both"
) + scale_color_manual(values = brewer.pal(6, "Paired")[3:6]) + scale_fill_manual(values = brewer.pal(6, "Paired")[3:6])
ggsave(
  "img/gsoap_gene_targets_pos_vs_neg_reg_mirnas_nolabel.pdf",
  width = 8,
  height = 6
)

#save
saveRDS(df, file = "data/go_terms_mir_targets_pos_and_neg.rds")

# --- MANUAL ANALYSIS OF GO TERM CLUSTERS IN PDF OUTPUT ####
df <- readRDS(file = "data/go_terms_mir_targets_pos_and_neg.rds")

# FURTHER ANALYSIS OF MANUAL CLUSTERS - COMMON GENES? ####
head(df)
#melt gene ensgs
df <- df %>% mutate(genes = strsplit(as.character(genes), "/")) %>% unnest(genes)
df %>% group_by(group) %>% dplyr::count(genes, sort = T) #top genes per group overall

#>reimport GO term clusters####
#paste from illustrator
#1: gene silencing by miRNA, mRNA processing, telomere, megakaryocytes, IL-7
r1 <- "GO:0035196: production of miRNAs involved in gene silencing by miRNAGO:1903313: positive regulation of mRNA metabolic processGO:0030219: megakaryocyte differentiationGO:0034249: negative regulation of cellular amide metabolic processGO:0060147: regulation of posttranscriptional gene silencingGO:0060968: regulation of gene silencingGO:0000291: nuclear−transcribed mRNA catabolic process, exonucleolyticGO:0010586: miRNA metabolic processGO:0000288: nuclear−transcribed mRNA catabolic process, deadenylation-dependent decay GO:0061014: positive regulation of mRNA catabolic processGO:0045814: negative regulation of gene expression  GO:0006334: nucleosome assembly GO:0032200: telomere organization GO:0000183: rDNA heterochromatin assembly GO:0051290: protein heterotetramerization GO:0035195: gene silencing by miRNA GO:0045652: regulation of megakaryocyte differentiationGO:0038111: interleukin−7−mediated signaling pathway "
p1 <- str_extract_all(r1, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#2: proliferation vs differentiation (de-repressed): cardiac muscle, neuronal (cortex), bone (BMP), epithelia; myelination, microtubules, NF-KB
r2 <- "GO:0045445: myoblast differentiationGO:0035914: skeletal muscle cell differentiationGO:0032409: regulation of transporter activityGO:0090090: negative regulation of canonical Wnt signaling pathwayGO:0043484: regulation of RNA splicingGO:0030279: negative regulation of ossificationGO:0010613: positive regulation of cardiac muscle hypertrophyGO:0050931: pigment cell differentiationGO:0060350: endochondral bone morphogenesisGO:0030155: regulation of cell adhesionGO:0042633: hair cycleGO:0055017: cardiac muscle tissue growthGO:0051149: positive regulation of muscle cell differentiationGO:0002062: chondrocyte differentiationGO:0022612: gland morphogenesisGO:0050680: negative regulation of epithelial cell proliferationGO:0008284: positive regulation of cell population proliferationGO:0019226: transmission of nerve impulseGO:0030510: regulation of BMP signaling pathwayGO:0006936: muscle contractionGO:0048873: homeostasis of number of cells within a tissueGO:0045596: negative regulation of cell differentiationGO:0032409: regulation of transporter activityGO:2000725: regulation of cardiac muscle cell differentiationGO:0010613: positive regulation of cardiac muscle hypertrophyGO:0048660: regulation of smooth muscle cell proliferationGO:0048646: anatomical structure formation involved in morphogenesisGO:0001508: action potentialGO:0055002: striated muscle cell developmentGO:0060038: cardiac muscle cell proliferationGO:0006367: transcription initiation from RNA polymerase II promoterGO:0010613: positive regulation of cardiac muscle hypertrophy"
p2 <- str_extract_all(r2, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#3 (de-repressed): B cells, T cells, cytokines, IFN-beta, defense response, LPS, NF-KB
r3 <- "GO:0002833: positive regulation of response to biotic stimulus
GO:2000108: positive regulation of leukocyte apoptotic processGO:0002363: alpha−beta T cell lineage commitmentGO:0032103: positive regulation of response to external stimulusGO:0070232: regulation of T cell apoptotic processGO:0046006: regulation of activated T cell proliferationGO:0000122: negative regulation of transcription by RNA polymerase IIGO:1901224: positive regulation of NIK/NF-kappaB signalingGO:2001141: regulation of RNA biosynthetic processGO:2000515: negative regulation of CD4-positive, alpha-beta T cell activation
GO:0032728: positive regulation of interferon-beta productionGO:0031349: positive regulation of defense responseGO:0032653: regulation of interleukin−10 productionGO:0031663: lipopolysaccharide−mediated signaling pathway
GO:0001959: regulation of cytokine−mediated signaling pathway
GO:0001819: positive regulation of cytokine productionGO:0009615: response to virusGO:0072538: T−helper 17 type immune responseGO:0043369: CD4−positive or CD8−positive alpha-beta T cell lineage commitmentGO:0032620: interleukin−17 productionGO:0001782: B cell homeostasisGO:0000122: negative regulation of transcription by RNA polymerase IIGO:0050853: B cell receptor signaling pathwayGO:0030888: regulation of B cell proliferationGO:0042093: T−helper cell differentiationGO:0030890: positive regulation of B cell proliferation"
p3 <- str_extract_all(r3, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#4 (de-repressed): response to stimulus: decreased oxygen levels, antibiotic, toxin, organophosphorus, glucocorticoid, nitrogen-containing, vitamin
r4 <- "GO:0001841: neural tube formationGO:0014009: glial cell proliferationGO:0050775: positive regulation of dendrite morphogenesisGO:0042552: myelinationGO:0042551: neuron maturationGO:0021545: cranial nerve developmentGO:1900449: regulation of glutamate receptor signaling pathwayGO:0048709: oligodendrocyte differentiationGO:0042490: mechanoreceptor differentiationGO:0021544: subpallium developmentGO:0007616: long−term memoryGO:0006998: nuclear envelope organizationGO:0032886: regulation of microtubule−based processGO:0007405: neuroblast proliferationGO:0001764: neuron migrationGO:0048839: inner ear developmentGO:0072384: organelle transport along microtubuleGO:0000186: activation of MAPKK activityGO:0021543: pallium developmentGO:0045773: positive regulation of axon extensionGO:0031646: positive regulation of nervous system processGO:0060977: coronary vasculature morphogenesisGO:0099518: vesicle cytoskeletal traffickingGO:1901880: negative regulation of protein depolymerizationGO:0036035: osteoclast developmentGO:0097435: supramolecular fiber organizationGO:0048813: dendrite morphogenesisGO:0060999: positive regulation of dendritic spine developmentGO:0031114: regulation of microtubule depolymerizationGO:0051491: positive regulation of filopodium assemblyGO:0048286: lung alveolus developmentGO:0009791: post−embryonic developmentGO:0035315: hair cell differentiationGO:0043547: positive regulation of GTPase activityGO:0001843: neural tube closureGO:0051056: regulation of small GTPase mediated signal transductionGO:0043123: positive regulation of I−kappaB kinase/NF-kappaB signalingGO:0002011: morphogenesis of an epithelial sheetGO:0002066: columnar/cuboidal epithelial cell developmentGO:0031098: stress−activated protein kinase signaling cascadeGO:0045197: establishment or maintenance of epithelial cell apical/basal polarityGO:0003401: axis elongationGO:0007266: Rho protein signal transductionGO:0030010: establishment of cell polarityGO:0120034: positive regulation of plasma membrane bounded cell projection assembly
GO:0051489: regulation of filopodium assembly "
p4 <- str_extract_all(r4, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#5: cellular location: ubiquitinylation, endosome/lysosome; cell cycle, PI3K, histone modification, NF-KB
r5 <- "GO:0006261: DNA−dependent DNA replicationGO:0046685: response to arsenic−containing substanceGO:0010390: histone monoubiquitinationGO:0031146: SCF−dependent proteasomal ubiquitin-dependent protein catabolic processGO:1903320: regulation of protein modification by small protein conjugation or removalGO:0000209: protein polyubiquitinationGO:0000731: DNA synthesis involved in DNA repairGO:0070646: protein modification by small protein removalGO:0031062: positive regulation of histone methylationGO:1902808: positive regulation of cell cycle G1/S phase transitionGO:0045023: G0 to G1 transitionGO:1901992: positive regulation of mitotic cell cycle phase transitionGO:0043124: negative regulation of I−kappaB kinase/NF-kappaB signalingGO:0008333: endosome to lysosome transportGO:0071260: cellular response to mechanical stimulusGO:0046839: phospholipid dephosphorylationGO:0001824: blastocyst developmentGO:2000679: positive regulation of transcription regulatory region DNA bindingGO:0034405: response to fluid shear stressGO:1904036: negative regulation of epithelial cell apoptotic process
GO:0045814: negative regulation of gene expression GO:0006644: phospholipid metabolic processGO:0014068: positive regulation of phosphatidylinositol 3-kinase signalingGO:0016567: protein ubiquitinationGO:0035306: positive regulation of dephosphorylationGO:0051443: positive regulation of ubiquitin−protein transferase activity"
p5 <- str_extract_all(r5, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#6 (de-repressed): response to stimulus: decreased oxygen levels, antibiotic, toxin, organophosphorus, glucocorticoid, nitrogen-containing, vitamin
r6 <- "GO:0036293: response to decreased oxygen levelsGO:1901699: cellular response to nitrogen compoundGO:0042493: response to drugGO:0007610: behaviorGO:0046683: response to organophosphorusGO:0033273: response to vitaminGO:0031668: cellular response to extracellular stimulusGO:0071453: cellular response to oxygen levelsGO:0051384: response to glucocorticoidGO:0010243: response to organonitrogen compoundGO:0097237: cellular response to toxic substanceGO:1901652: response to peptideGO:0071236: cellular response to antibioticGO:0071417: cellular response to organonitrogen compoundGO:0036293: response to decreased oxygen levelsGO:0031669: cellular response to nutrient levels"
p6 <- str_extract_all(r6, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#7 (repressed): IL-1 and IL-6 production, humoral antimicrobial immune response, prostaglandins
r7 <- "GO:0032755: positive regulation of interleukin−6 productionGO:0002690: positive regulation of leukocyte chemotaxisGO:0034694: response to prostaglandinGO:0060760: positive regulation of response to cytokine stimulusGO:0002720: positive regulation of cytokine production involved in immune responseGO:0032732: positive regulation of interleukin−1 productionGO:0031640: killing of cells of other organism GO:0006693: prostaglandin metabolic process GO:0006636: unsaturated fatty acid biosynthetic process GO:0050832: defense response to fungus GO:0006805: xenobiotic metabolic process GO:0046456: icosanoid biosynthetic process GO:0019369: arachidonic acid metabolic process GO:0061844: antimicrobial humoral immune response  mediated by antimicrobial peptideGO:0050701: interleukin−1 secretion GO:0006959: humoral immune response "
p7 <- str_extract_all(r7, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#8: ion transport and electrical coupling, sodium, potassium, calcium, blood circulation and cardiac action (mixed), pain, insulin (mixed)
r8 <- "GO:0098901: regulation of cardiac muscle cell actionGO:0050796: regulation of insulin secretionGO:0090279: regulation of calcium ion importGO:0010644: cell communication by electrical couplingGO:0010765: positive regulation of sodium ion transportGO:0043268: positive regulation of potassium ion transport GO:0035725: sodium ion transmembrane transportGO:0050796: regulation of insulin secretionGO:0070588: calcium ion transmembrane transportGO:0008277: regulation of G protein−coupled receptorGO:0015837: amine transportGO:1903523: negative regulation of blood circulationGO:0051930: regulation of sensory perception of painGO:0007204: positive regulation of cytosolic calciumGO:0060314: regulation of ryanodine−sensitive calcium-release channel activityGO:1990573: potassium ion import across plasma membrane
GO:0071277: cellular response to calcium ion "
p8 <- str_extract_all(r8, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#9: synaptic processes, transport, motility, autophagy, lipoprotein and phospholipid
r9 <- "GO:0046794: transport of virusGO:0048278: vesicle dockingGO:0006892: post−Golgi vesicle−mediated transportGO:0046794: transport of virusGO:0098693: regulation of synaptic vesicle cycleGO:0010508: positive regulation of autophagyGO:0090150: establishment of protein localization to membraneGO:0035418: protein localization to synapseGO:0042158: lipoprotein biosynthetic processGO:0061025: membrane fusionGO:0045332: phospholipid translocationGO:0060285: cilium−dependent cell motility GO:0048168: regulation of neuronal synaptic plasticity"
p9 <- str_extract_all(r9, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#10: lipid and glucose metabolism, cardiac muscle apoptosis, steroid hormone regulation (stress), EGFR
r10 <- "GO:0019915: lipid storageGO:0016052: carbohydrate catabolic processGO:0050873: brown fat cell differentiationGO:0060669: embryonic placenta morphogenesisGO:0046885: regulation of hormone biosynthetic processGO:0140353: lipid export from cellGO:0071402: cellular response to lipoprotein particle stimulus
GO:0007176: regulation of epidermal growth factor-activated receptor activityGO:0055088: lipid homeostasisGO:0097066: response to thyroid hormoneGO:0008209: androgen metabolic processGO:0043620: regulation of DNA−templated transcription in response to stressGO:0015850: organic hydroxy compound transportGO:0045821: positive regulation of glycolytic processGO:0006641: triglyceride metabolic processGO:0055081: anion homeostasisGO:0045923: positive regulation of fatty acid metabolic processGO:0006493: protein O−linked glycosylationGO:0046058: cAMP metabolic processGO:0002673: regulation of acute inflammatory responseGO:0050994: regulation of lipid catabolic processGO:1905954: positive regulation of lipid localizationGO:0010659: cardiac muscle cell apoptotic processGO:0010656: negative regulation of muscle cell apoptotic processGO:0032007: negative regulation of TOR signalingGO:1905039: carboxylic acid transmembrane transportGO:0009100: glycoprotein metabolic processGO:0045834: positive regulation of lipid metabolic processGO:0019217: regulation of fatty acid metabolic processGO:0010828: positive regulation of glucose transmembrane transportGO:0034198: cellular response to amino acid starvationGO:0006654: phosphatidic acid biosynthetic processGO:0032196: transposition GO:0008207: C21−steroid hormone metabolic process "
p10 <- str_extract_all(r10, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#11: cell cycle checkpoint, DNA damage control
r11 <- "GO:0072395: signal transduction involved in cell cycle checkpointGO:0007093: mitotic cell cycle checkpointGO:0000077: DNA damage checkpointGO:0007050: cell cycle arrest"
p11 <- str_extract_all(r11, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#12: cholesterol synthesis, thermogenesis
r12 <- "GO:0120161: regulation of cold−induced thermogenesisGO:0006084: acetyl−CoA metabolic processGO:0090181: regulation of cholesterol metabolic processGO:0010675: regulation of cellular carbohydrate metabolic processGO:1902930: regulation of alcohol biosynthetic processGO:0016125: sterol metabolic process "
p12 <- str_extract_all(r12, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

#13: chromatin organisation, gene silencing (repressed), gene expression (de-repressed)
r13 <- "GO:0045815: positive regulation of gene expression  GO:1905268: negative regulation of chromatin organizationGO:0031935: regulation of chromatin silencing GO:0060969: negative regulation of gene silencing "
p13 <- str_extract_all(r13, "GO:(.*?):") %>% unlist() %>% substr(1, 10)

term_clusters <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13)
names(term_clusters) <- paste0("cluster", 1:13)
term_genes <- vector(mode = "list", length = length(term_clusters))
names(term_genes) <- names(term_clusters)

#>foreach cluster analysis####
for (i in 1:length(term_genes)) {
  cl <- term_clusters[[i]]
  df_t <- df[df$GO.ID %in% cl,]
  top_genes <- df_t %>% group_by(group) %>% dplyr::count(genes, sort = T)
  top_genes$symbol <- ensg2symbol$hgnc_symbol[match(top_genes$genes, ensg2symbol$ensembl_gene_id)]
  top_genes$perc <- top_genes$n / length(cl)
  term_genes[[i]] <- top_genes
}
term_genes
top_gene_df <- plyr::ldply(term_genes)
top_gene_df <- top_gene_df[order(top_gene_df$n, decreasing = T),]
nrow(top_gene_df)
head(top_gene_df, 100)
#account for different cluster sizes
#percent?
hist(top_gene_df$perc)
quantile(top_gene_df$perc, seq(0,1,.01))

table(top_gene_df$.id)

#>> per cluster to address size differences####
nrow(df)
cl_list <- lapply(df$GO.ID, function(id) which(unlist(lapply(term_clusters, function(x) any(grepl(id, x))))))
df$cluster <- lapply(cl_list, paste, collapse = ", ") %>% unlist()
df <- df %>% mutate(cluster = strsplit(as.character(cluster), ", ")) %>% unnest(cluster)
df$cluster <- paste0("cluster", df$cluster)
df <- df %>% mutate(ensg = genes) %>% dplyr::select(GO.ID, ensg, cluster) %>% distinct()
df$symbol <- ensg2symbol$hgnc_symbol[match(df$ensg, ensg2symbol$ensembl_gene_id)]

#fisher exact test
df_top <- df %>% group_by(cluster) %>% dplyr::count(ensg, sort = T) %>% ungroup()
df_top$pval <- NA

#loop over
for (i in 1:nrow(df_top)) {
  m <- matrix(c(nrow(df[df$ensg == df_top$ensg[i] & df$cluster == df_top$cluster[i],]),
              nrow(df[df$ensg == df_top$ensg[i] & df$cluster != df_top$cluster[i],]),
              nrow(df[df$ensg != df_top$ensg[i] & df$cluster == df_top$cluster[i],]),
              nrow(df[df$ensg != df_top$ensg[i] & df$cluster != df_top$cluster[i],])),
              nrow = 2, ncol = 2)
  t <- fisher.test(m, alternative = "greater")
  df_top$pval[i] <- t$p.value
}
df_top$symbol <- ensg2symbol$hgnc_symbol[match(df_top$ensg, ensg2symbol$ensembl_gene_id)]
df_top <- df_top[order(df_top$n, decreasing = T), ]

# reorder columns
df_top <- df_top %>% dplyr::select(ensg, symbol, everything())

#adjust pval
df_top$padj <- p.adjust(df_top$pval, method = "BH")
df_top$de_padj <- de_genes$padj[match(df_top$ensg, de_genes$ensg)]
df_top$de_baseMean <- de_genes$baseMean[match(df_top$ensg, de_genes$ensg)]
df_top$de_logFC <- de_genes$log2FoldChange[match(df_top$ensg, de_genes$ensg)]
df_top$de_countChange <- de_genes$countChange[match(df_top$ensg, de_genes$ensg)]

#>>>write out####
write.csv(df_top, file = "out/mir_targets_go_cluster_enrichment.csv", quote = F, row.names = F)

mean(df_top$padj<.05)
mean(df_top$padj<.001)

df_top_05 <- data.frame(df_top[df_top$padj<.05, ])
nrow(df_top_05)

df_top_05$cluster <- gsub("cluster", "", df_top_05$cluster)
cl_count <- dplyr::count(df_top_05, cluster, wt = n(), sort = T) #wt = n() or else it weighs by n column
df_top_05$cluster <- factor(df_top_05$cluster, levels = rev(cl_count$cluster))

ggplot(df_top_05, aes(cluster)) + geom_bar() + coord_flip() +
  theme_minimal()
ggsave("img/mir_target_go_enriched_genes_per_cluster.pdf", width = 5, height = 7)

df_top_cl6 <- df_top_05 %>% filter(cluster == "6")
df_top_cl8 <- df_top_05 %>% filter(cluster == "8")

#check DAVID
write.table(df_top_05$ensg, file = "out/mir_targets_go_all_below_05.txt", sep = "\n", row.names = F, quote = F)
write.table(df_top_cl6$ensg, file = "out/mir_targets_go_cluster6.txt", sep = "\n", row.names = F, quote = F)
write.table(df_top_cl8$ensg, file = "out/mir_targets_go_cluster8.txt", sep = "\n", row.names = F, quote = F)
