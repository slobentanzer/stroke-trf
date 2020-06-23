# NETWORK GENERATION WAS PERFORMED USING MIRNEO (INHOUSE DATABASE, see Lobentanzer et al. 2019, Cell Reports) #
# The CD14-specific network of differentially expressed transcription factors (Figure 5B) #
# was exported as .gexf and spatialised via gephi 0.9. (File: "out/DE_TFs_miRs_tRFs_CD14.gexf") #
# The same network was used to calculate miRNA/tRF targeting ratios (Figure 5 C&D). #

setwd("~/GitHub/stroke-trf")
rm(list = ls())

library(ggplot2)

lnames <- load("data/CD14_TF_miR_tRF_network.RData")
lnames

#BAR GRAPH OF TFs TARGETING MIR vs TRF FRACTIONS####
find.fraction <- function(id){
  if(id %in% edges$target) {
    temp <- edges[edges$target == id, ]
    temp$biotype <- nodes$biotype[match(temp$source, nodes$id)]
    return(mean(temp$biotype == "tRF"))
  } else {
    return(2)
  }
}

nodes$tRF_fraction <- unlist(lapply(nodes$id, find.fraction))
nodes[nodes$biotype == "TF_top18",]

fnodes <- nodes[!is.na(nodes$tRF_fraction),]
fnodes <- fnodes[grep("TF", fnodes$biotype),]

fnodes <- fnodes[!is.na(fnodes$cc),]
fnodes <- fnodes[order(abs(fnodes$cc), decreasing = T),]
fnodes$tRF_fraction <- fnodes$tRF_fraction - .5

fnodes$log_cc <- log2(abs(fnodes$cc) + 1)
nrow(fnodes)

#>cholinergic TFs?####
cholinergic_tfs <- readRDS("data/cholinergic_tfs.rds")
cholinergic_genes <- readRDS("data/cholinergic_genes.rds")
cholinergic_tfs <- cholinergic_tfs[cholinergic_tfs$tissue == "CD14POS_MONOCYTES",]
cholinergic_tfs <- cholinergic_tfs[cholinergic_tfs$g.name %in% cholinergic_genes$gene_symbol[cholinergic_genes$group %in% c("core", "receptor")],]
unique(cholinergic_tfs$tf.name[cholinergic_tfs$g.name == "SLC18A3"])
fnodes$label %in% cholinergic_tfs$tf.name
fnodes$cholinergic <- fnodes$label %in% cholinergic_tfs$tf.name
table(fnodes$cholinergic)

fnodes <- fnodes[order(fnodes$log_cc),]
fnodes <- fnodes[order(fnodes$tRF_fraction),]
fnodes$label <- factor(fnodes$label, levels = unique(fnodes$label))

#>>display all targeted tfs####
nrow(fnodes[fnodes$tRF_fraction<1.5,])
ggplot(fnodes[fnodes$tRF_fraction<1.5,], aes(label, tRF_fraction, fill = tRF_fraction)) + 
  geom_bar(stat = "identity") + 
  scale_fill_viridis_c() +
  coord_flip() + #ggtitle(paste0("tRF fraction of top TFs")) +
  ylab("tRF fraction - 0.5") + xlab("TF name") +
  theme_minimal()
ggsave("img/TF_trf_fraction_all_bar.svg", width = 8, height = 5)

#handle non-targeted TFs
fnodes$targeted <- fnodes$tRF_fraction<1
table(fnodes$targeted)
table(fnodes$tRF_fraction[fnodes$targeted] == .5) #only trfs
table(fnodes$tRF_fraction[fnodes$targeted] == -.5) #only mirs

fnodes$tRF_fraction[!fnodes$targeted] <- fnodes$tRF_fraction[!fnodes$targeted] - 1.5
fnodes <- fnodes[order(fnodes$targeted, fnodes$tRF_fraction),]
fnodes$label <- factor(fnodes$label, levels = unique(fnodes$label))

fnodes1 <- fnodes[fnodes$padj < .1,]

ggplot(fnodes1, aes(label, tRF_fraction, col = cc>0, size = log_cc)) + 
  geom_point(stat = "identity") + 
  scale_color_discrete(name = "up-regulated") +
  scale_size(range = c(1,7), name = "count change", breaks = fivenum(fnodes1$log_cc[1:20]), 
             limits = c(min(fnodes1$log_cc[1:20]), max(fnodes1$log_cc[1:20])),
             labels = round(2^fivenum(fnodes1$log_cc[1:20]), 0)) +
  geom_text(aes(fnodes1$label[which(fnodes1$label[1:20] %in% cholinergic_tfs$tf.name)], 
                fnodes1$tRF_fraction[which(fnodes1$label[1:20] %in% cholinergic_tfs$tf.name)]), 
            label = "c", position = "identity", 
            data = fnodes1[which(fnodes1$label[1:20] %in% cholinergic_tfs$tf.name),],
            col = "black") + 
  coord_flip() + #ggtitle(paste0("tRF fraction of top TFs")) +
  ylab("tRF fraction - 0.5") + xlab("TF name") +
  scale_y_continuous(limits = c(-.6,.6), breaks = c(-.5,-.25,0,.25,.5)) +
  theme_minimal()
ggsave("img/TF_tRF_fraction_top26_alpha0.1.svg", width = 4, height = 8)


fnodes0 <- fnodes[fnodes$padj < .05,]
nrow(fnodes0)

ggplot(fnodes0, aes(label, tRF_fraction, col = cc>0, size = log_cc)) + 
  geom_point(stat = "identity") + 
  scale_color_discrete(name = "up-regulated") +
  scale_size(range = c(3,10), name = "count change", breaks = fivenum(fnodes0$log_cc[1:20]), 
             limits = c(min(fnodes0$log_cc[1:20]), max(fnodes0$log_cc[1:20])),
             labels = round(2^fivenum(fnodes0$log_cc[1:20]), 0)) +
  geom_text(aes(fnodes0$label[which(fnodes0$label %in% cholinergic_tfs$tf.name)], 
                fnodes0$tRF_fraction[which(fnodes0$label %in% cholinergic_tfs$tf.name)]), 
            label = "c", position = "identity", 
            data = fnodes0[which(fnodes0$label %in% cholinergic_tfs$tf.name),],
            col = "black") + 
  coord_flip() + #ggtitle(paste0("tRF fraction of top TFs")) +
  ylab("tRF fraction - 0.5") + xlab("TF name") +
  scale_y_continuous(limits = c(-.6,.6), breaks = c(-.5,-.25,0,.25,.5)) +
  theme_minimal()
ggsave("img/TF_tRF_fraction_top18_alpha0.05.svg", width = 4, height = 8)
