setwd("~/GitHub/stroke-trf")
rm(list = ls())

### code exported from BioRad CFX Maestro Software version 4.1.24.33.1219

# experiments with ssRNA tRF-22-WE8SPOX52 mimics
# quantification of its target Zbp1 24h after transfection of RAW 264.7 cells (Figure 6F)

#check and install required packages
pacman::p_load(car, lsmeans)

ShapiroTest <- function(x){
  z <- length(x)
  if(z > 2 && z < 5001){
    y <- shapiro.test(x)
    return (c(w=unname(y$statistic), p=y$p.value, n=z))
  }
  return (c(w=NA, p=NA, n=z))
}

OneWayAnova <- function(data, dv, iv, conf = 0.95, toJSON = FALSE){
 
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

#Zbp-1
expzbp <- c(-0.662567617782518,-0.4490718770484,-1.10776210348588,-0.119619190804873,0.135262015952943,-0.015642825148067)
bgzbp <- c("mimics","mimics","mimics","NC","NC","NC")
dfzbp <- data.frame(expzbp, bgzbp)
OneWayAnova(data = dfzbp, dv = "expzbp", iv = "bgzbp", conf = 0.95)

# quantification of the tRF-22-WE8SPOX52 and its parental tRNA to confirm overexpression
# (Fig S7)

# experiment 1 - used for RNA sequencing

#tRF-22-WE8SPOX52
exptrf.22 <- c(4.83064593635308,4.48306365898003,4.3326768556368,-0.573442782462649,0.368738330917871,0.204704451544778)
bgtrf.22 <- c("mimics 22","mimics 22","mimics 22","NC","NC","NC")
dftrf.22 <- data.frame(exptrf.22, bgtrf.22)
OneWayAnova(data = dftrf.22, dv = "exptrf.22", iv = "bgtrf.22", conf = 0.95)

# experiment 2 - used for RT-qPCR validations

#tRF-22-WE8SPOX52
exptrf22 <- c(1.52144683930338,2.7008483312321,2.03092891397204,0.0502232405876626,-0.234576421318813,0.184353180731142)
bgtrf22 <- c("mimics","mimics","mimics","NC","NC","NC")
dftrf22 <- data.frame(exptrf22, bgtrf22)
OneWayAnova(data = dftrf22, dv = "exptrf22", iv = "bgtrf22", conf = 0.95)