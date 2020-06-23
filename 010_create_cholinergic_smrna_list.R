setwd("~/GitHub/stroke-trf")
rm(list = ls())

#CHOLINO-MIRS AND -TRFS####
#load graph output
cholinergic_trfs <- readRDS("data/cholinergic_trfs.rds")
cholinergic_mirs <- readRDS("data/cholinergic_mirs.rds")

#find smRNAs targeting multiple cholinergic transcripts
####mir
mir_count <- dplyr::count(cholinergic_mirs[cholinergic_mirs$r.rating>5 & !cholinergic_mirs$group %in% c("extended", "nt", "nt_receptor", "circadian"),], m.name, sort = T)
write.csv(mir_count, "out/cholinomiR_counts.csv", row.names = F, quote = F)
plot(ecdf(mir_count$n))
abline(h = .8, col = "red")
quantile(mir_count$n)
quantile(mir_count$n, .8)
pdf("img/cholinomir_ecdf.pdf", height = 8, width = 8)
plot(ecdf(mir_count$n))
abline(h = .8, col = "red")
dev.off()

nrow(mir_count[mir_count$n>4,])
cholinomirs <- mir_count$m.name[mir_count$n>4]

####trf
nrow(unique(cholinergic_trfs[, c("MINTplate", "ensg")]))
trf_count <- dplyr::count(unique(cholinergic_trfs[!cholinergic_trfs$group %in% c("extended", "nt", "nt_receptor", "circadian"), c("MINTplate", "ensg")]), MINTplate, sort = T)
write.csv(trf_count, "out/cholinotrf_counts.csv", row.names = F, quote = F)
plot(ecdf(trf_count$n))
quantile(trf_count$n)
quantile(trf_count$n, .8)
abline(h = .8, col = "red")
pdf("img/cholinotrf_ecdf.pdf", height = 8, width = 8)
plot(ecdf(trf_count$n))
abline(h = .8, col = "red")
dev.off()

nrow(trf_count[trf_count$n>4,])
cholinotrfs <- trf_count$MINTplate[trf_count$n>4]

save(cholinomirs, cholinotrfs, file = "data/cholinergic_associated_smrnas.RData")
