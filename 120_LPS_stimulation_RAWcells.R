setwd("~/GitHub/stroke-trf")
rm(list = ls())

### code exported from BioRad CFX Maestro Software version 4.1.24.33.1219

#LPS stimualation experiments in murine RAW 264.7 cells: 

#quantification of the top 6 stroke-perturbed tRFs 18h after LPS (Figure 6 A&B)

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

#tRF-18-8R6546D2
exp188R5 <- c(0.288575142144076,-0.358460824534921,0.0698856823908526,1.09540500025325,0.882883015905381,0.902190865192479,0.193789012287732,-0.182041720423483)
bg188R5 <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
df188R5 <- data.frame(exp188R5, bg188R5)
OneWayAnova(data = df188R5, dv = "exp188R5", iv = "bg188R5", conf = 0.95)

#tRF-18-HR0VX6D2
exp18HRO <- c(0.335680446656038,-0.25030761113222,-0.085372835523808,0.788779294180509,0.604536511922763,0.448042925867723,0.192338803278025,0.128334161976968)
bg18HRO <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
df18HRO <- data.frame(exp18HRO, bg18HRO)
OneWayAnova(data = df18HRO, dv = "exp18HRO", iv = "bg18HRO", conf = 0.95)

#tRF-22-8EKSP1852
exp228EK <- c(0.242059958306008,-0.28349365444649,0.0414336961404818,1.0047225767457,0.842697672121431,0.726693207454039,0.227358747688017,0.272596414646566)
bg228EK <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
df228EK <- data.frame(exp228EK, bg228EK)
OneWayAnova(data = df228EK, dv = "exp228EK", iv = "bg228EK", conf = 0.95)

#tRF-22-WE8SPOX52
exp22SPOX <- c(0.388071726743114,-0.457968238441152,0.0698965116980385,1.19493036411597,0.964059918030767,0.964178592617651,0.503762458591799,0.674944362673868)
bg22SPOX <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
df22SPOX <- data.frame(exp22SPOX, bg22SPOX)
OneWayAnova(data = df22SPOX, dv = "exp22SPOX", iv = "bg22SPOX", conf = 0.95)

#tRF-22-WEKSPM852
exp22WEK <- c(0.319957851400137,-0.288197040495897,-0.0317608109042294,1.09528018891663,0.875065607005322,0.916486846487752,0.373865405601902,0.373893531228445)
bg22WEK <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
df22WEK <- data.frame(exp22WEK, bg22WEK)
OneWayAnova(data = df22WEK, dv = "exp22WEK", iv = "bg22WEK", conf = 0.95)

#tRF-18-8R6Q46D2
exp188RQ <- c(0.123617156530237,-0.289937754813293,0.166320598283059,0.917947842722198,0.956859817626885,0.754669757519741,0.170709570330837,0.350832620055394)
bg188RQ <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
df188RQ <- data.frame(exp188RQ, bg188RQ)
OneWayAnova(data = df188RQ, dv = "exp188RQ", iv = "bg188RQ", conf = 0.95)

# expression of inflammatory signaling moleculaes 18h after LPS (Figure S6)

#Cd14
expcd14 <- c(0.196352382546687,0.0718030386920763,-0.268155421238754,2.51210719548195,2.63938558365553,2.45938516023407,1.96717227452414,1.52706009222983)
bgcd14 <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
dfcd14 <- data.frame(expcd14, bgcd14)
OneWayAnova(data = dfcd14, dv = "expcd14", iv = "bgcd14", conf = 0.95)

#Il10
expil10 <- c(-0.122741617151403,0.314046327624856,-0.191304710473448,4.92380400995986,4.84124763781024,3.74896005291058,1.72633189806664,0.235752976166948)
bgil10 <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
dfil10 <- data.frame(expil10, bgil10)
OneWayAnova(data = dfil10, dv = "expil10", iv = "bgil10", conf = 0.95)

#Stat1
expstat1 <- c(0.256710047962214,0.226356956442347,-0.483067004404558,3.23808915799921,2.59248875675332,2.64308775253871,0.492318785326808,0.220190982635925)
bgstat1 <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
dfstat1 <- data.frame(expstat1, bgstat1)
OneWayAnova(data = dfstat1, dv = "expstat1", iv = "bgstat1", conf = 0.95)

#Tnfa
exptnfa <- c(-0.00792193800583467,0.00090291132640115,0.00701902667944267,3.48841581336172,3.62606392341475,3.0377343237112,2.81516026652955,1.87148056516409)
bgtnfa <- c("nonstimulated","nonstimulated","nonstimulated","stimulated","stimulated","stimulated","dexamethasone","dexamethasone")
dftnfa <- data.frame(exptnfa, bgtnfa)
OneWayAnova(data = dftnfa, dv = "exptnfa", iv = "bgtnfa", conf = 0.95)