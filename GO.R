#TARGETS GO ANALYSIS####
## >>unbiased GO analysis - all results####
library(topGO)

get.GO <- function(scores, ensgs, quant = .9, pval = F, p.threshold = .05, top = 2000,
                   topGenes = function(allScore){return(allScore >= x)}){
  geneID2GO <- readRDS(file = "data/geneID2GO_all.rds")
  
  #ontology
  ontology <- "BP"
  
  allGenes <- as.numeric(scores[1:top])
  names(allGenes) <- ensgs[1:top]
  
  hist(allGenes)
  if(!pval) {
    x <- quantile(allGenes, seq(0, 1, quant))
    topGenes <- topGenes
  } else {
    topGenes <- function(allScore){return(allScore < p.threshold)}
  }
  
  
  #prune GO terms
  nodeSize <- 10
  
  GOdata <- new(
    "topGOdata",
    # description = description,
    ontology = ontology,
    allGenes = allGenes,
    geneSel = topGenes,
    annot = annFUN.gene2GO,
    nodeSize = nodeSize,
    gene2GO = geneID2GO
  )
  
  ## >fisher test ####
  test.stat <-
    new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.stat)
  resultFisher
  
  # >weight test ####
  test.stat <-
    new(
      "weightCount",
      testStatistic = GOFisherTest,
      name = "Fisher test",
      sigRatio = "ratio"
    )
  resultWeight <- getSigGroups(GOdata, test.stat)
  resultWeight
  
  ## >total ####
  allRes <- GenTable(
    GOdata,
    classic = resultFisher,
    weight = resultWeight,
    orderBy = "weight",
    ranksOf = "classic",
    topNodes = 100
  )
  (res <- allRes[allRes$weight < .05, ])
  
  return(list(result = res, GOdata = GOdata, fisher = resultFisher, weight = resultWeight))
}