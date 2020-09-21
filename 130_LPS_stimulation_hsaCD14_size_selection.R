setwd("~/GitHub/stroke-trf")
rm(list = ls())

### code exported from BioRad CFX Maestro Software version 4.1.24.33.1219

#LPS stimualation experiments in human CD14+ monocytes: 

#quantification of the top 6 stroke-perturbed tRFs 6h after LPS

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
exp188R5 <- c(-0.209549671544299,0.152721510035469,-0.0149301359625476,0.241162226021677,-0.35011229672649,0.180708368176186,-0.0673043020611762,0.156420554881256,0.221953334692266,0.207029301099883,0.485055343019168,0.498817273699633,-0.0037094409739461,0.508972007729486,0.421914881783653,-0.0490934284096464,0.403330111495267,0.701547783481922,-0.298595468675932,-0.342583180459298,0.0678597049991525,-0.0209142061559593,-0.294729713432565,0.216680452636121)
bg188R5 <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df188R5 <- data.frame(exp188R5, bg188R5)
OneWayAnova(data = df188R5, dv = "exp188R5", iv = "bg188R5", conf = 0.95)

#tRF-18-8R6Q46D2
exp188RQ <- c(-0.266208550616401,0.0127434518423964,-0.103163103491951,0.454083083119338,-0.247450311710892,0.149995430857514,-0.132775132607444,0.112078715283884,0.185066752094695,0.128804356725379,0.413152114506465,0.531373444805194,0.103900836816809,0.427991951860015,0.254656608165881,0.240834004329651,0.532059253423233,0.604556986226713,-0.19931410816327,-0.409059551770005,0.0152508462679832,0.107557457859629,-0.0189059156173367,0.336497855101049)
bg188RQ <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df188RQ <- data.frame(exp188RQ, bg188RQ)
OneWayAnova(data = df188RQ, dv = "exp188RQ", iv = "bg188RQ", conf = 0.95)

#tRF-18-HR0VX6D2
exp18HRO <- c(-0.575797025374705,-0.332299770083441,-0.578675191839089,0.0471780040815037,1.27423184687724,0.165362136338478,-0.621205870758884,-0.845324421782921,-0.173792094764543,-0.224437921104158,0.435718316238896,0.366235577927093,-0.111698576933193,-0.462225820261622,-0.530699907023291,-0.226433785823724,0.340209418461926,0.262869484748412,-0.270019440314476,-0.853737476887441,-0.315972232679022,-0.189284016363741,0.143948464141894,0.292271455259543)
bg18HRO <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df18HRO <- data.frame(exp18HRO, bg18HRO)
OneWayAnova(data = df18HRO, dv = "exp18HRO", iv = "bg18HRO", conf = 0.95)

#tRF-228EKSP1852
exp228EK <- c(-0.1086060254259,-0.0829055696082327,0.146488581050068,0.290404715259575,-0.206045417628289,-0.0393362836472146,0.0863477833480265,0.20553150017705,0.415735401414231,0.102768890401979,0.628909968703001,0.484491439590535,-0.0429326250452853,0.313691329335885,0.290487903947483,-0.111641373958184,0.554737457724999,0.359021512998984,-0.112944811285035,-0.174567042614134,0.262205539565018,0.0618845312582027,0.0408143825889001,0.131252369466052)
bg228EK <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df228EK <- data.frame(exp228EK, bg228EK)
OneWayAnova(data = df228EK, dv = "exp228EK", iv = "bg228EK", conf = 0.95)

#tRF-22-WE8SPOX52
exp22SPOX <- c(-0.440396118773613,-0.158960987595915,-0.20250570923973,0.0379371455219245,0.770837361507894,-0.00691169142056396,-0.20178956838066,0.071397712222671,0.337975975045986,-0.0969560529893662,0.408492067527918,0.361895662828285,-0.0940666828567985,0.233417843476669,0.102551230015702,-0.242401036854831,0.304092488348983,0.390365574811733,-0.225535314953982,-0.348834262859081,0.126059649504734,-0.0996295256556133,0.06374268445275,0.345027649257672)
bg22SPOX <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df22SPOX <- data.frame(exp22SPOX, bg22SPOX)
OneWayAnova(data = df22SPOX, dv = "exp22SPOX", iv = "bg22SPOX", conf = 0.95)

#tRF-22-WEKSPM852
exp22WEK <- c(-0.365932301823307,-0.174942865873877,-0.134593884048926,0.359641792651866,-0.00105145118753564,0.316878710281776,-0.262287723315717,-0.14655511646029,0.0698156962225909,0.0462050148685691,0.466188672988823,0.504992753252459,-0.0193717467672273,0.045467505941975,-0.00281524018639336,-0.000274545478026504,0.306690505136059,0.512175252281272,-0.295105445598741,-0.523753379585274,-0.17365621106306,-0.0310412608832084,-0.137802104490442,0.263697596491175)
bg22WEK <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df22WEK <- data.frame(exp22WEK, bg22WEK)
OneWayAnova(data = df22WEK, dv = "exp22WEK", iv = "bg22WEK", conf = 0.95)

#quantification of the top 6 stroke-perturbed tRFs 12h after LPS

#tRF-18-8R6546D2
exp188R5 <- c(-0.0829053627293704,-0.0626407560395386,-0.0926314764596529,-0.012984832200557,0.0779896121439946,0.173172815285111,0.336852792764994,0.313494900155581,0.402241258098163,0.542250038801072,0.375509103724877,0.446767280639094,0.321641099990943,0.627210428793864,0.725079399299629,0.71181249161836,0.44946482311023,0.479856287412044,0.0025561210350773,0.0829321110234762,-0.205382372493354,0.149951787029192,0.215564545874413,0.092527652328327)
bg188R5 <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df188R5 <- data.frame(exp188R5, bg188R5)
OneWayAnova(data = df188R5, dv = "exp188R5", iv = "bg188R5", conf = 0.95)

#tRF-18-8R6Q46D2
exp188RQ <- c(-0.174636647145512,-0.0789459937060747,-0.14348583829389,0.109267980344805,0.150821837738121,0.136978661062573,0.247787633805654,0.0989070896823403,0.0966091303496892,0.571054895400163,0.16511234376464,0.377987512608691,0.209981287436969,0.529182285274727,0.588204857449922,0.662660345655254,0.351022461863924,0.447953736551007,-0.128969809928456,0.0613244542780382,-0.314332571134958,0.342641606995389,0.191467861045842,0.0727020849266237)
bg188RQ <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df188RQ <- data.frame(exp188RQ, bg188RQ)
OneWayAnova(data = df188RQ, dv = "exp188RQ", iv = "bg188RQ", conf = 0.95)

#tRF-18-HR0VX6D2
exp18HRO <- c(0.0746742368531556,-0.366681872164442,-0.397031895262723,0.18587049356154,0.0961519734713254,0.407017063541138,0.158519091466388,-0.379243595514295,-0.377443479637845,0.275547811367997,0.0653041259900427,0.328174747132589,0.184156456007939,0.07450434121119,-0.13895124490451,0.79129201353852,0.28869173913229,0.396906098433238,0.22088485647884,0.000283630002739351,-0.473338426289093,0.472213466575553,0.208853882640609,0.228776988123789)
bg18HRO <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df18HRO <- data.frame(exp18HRO, bg18HRO)
OneWayAnova(data = df18HRO, dv = "exp18HRO", iv = "bg18HRO", conf = 0.95)

#tRF-228EKSP1852
exp228EK <- c(-0.0995852469041518,-0.0534874331022834,-0.0702116643627005,-0.0140640822984059,0.110971657612815,0.126376769054728,0.404521365798114,0.258185525647095,0.279923519612746,0.721211048823555,0.318475255121531,0.400564404511979,0.303354383113328,0.666410441897751,0.710430522409146,0.746032345232479,0.360608543318613,0.399360719826625,0.167095072545099,-0.0630902339714705,-0.0445265267830379,0.187907747820315,0.30581590939653,0.0472430104207806)
bg228EK <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df228EK <- data.frame(exp228EK, bg228EK)
OneWayAnova(data = df228EK, dv = "exp228EK", iv = "bg228EK", conf = 0.95)

#tRF-22-WE8SPOX52
exp22SPOX <- c(-0.100331101414476,0.00206260937552816,-0.142697497960517,0.110430287860508,-0.133374542831142,0.263910244970108,0.205518415964161,0.198909916290345,0.281372916851661,0.515480248713736,0.128223025523141,0.329073238557694,0.236508483697841,0.623091161367229,0.570446842845994,0.671763537971159,0.394210649996396,0.444898855974408,0.064406137639608,0.150170901426877,-0.0464409092086553,0.336542699635327,0.366611017870279,0.0128848413753583)
bg22SPOX <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df22SPOX <- data.frame(exp22SPOX, bg22SPOX)
OneWayAnova(data = df22SPOX, dv = "exp22SPOX", iv = "bg22SPOX", conf = 0.95)

#tRF-22-WEKSPM852
exp22WEK <- c(-0.114388972804065,-0.0886980880425625,-0.198765013230409,0.0446902012340169,0.176350367031368,0.180811505811652,0.279334741582371,0.100054740331954,0.17230555318457,0.583985913398813,0.365690337838821,0.500119357744939,0.227215577718018,0.545055354018872,0.575744219150702,0.824211746421433,0.561800858996439,0.536611783554782,-0.026218903557215,-0.00673178113931357,-0.381743090894581,0.29802122511477,0.24749755376852,0.0166646606512392)
bg22WEK <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df22WEK <- data.frame(exp22WEK, bg22WEK)
OneWayAnova(data = df22WEK, dv = "exp22WEK", iv = "bg22WEK", conf = 0.95)

#quantification of the top 6 stroke-perturbed tRFs 18h after LPS

#tRF-18-8R6546D2
exp188R5 <- c(-0.352503934890473,0.249169880145363,-0.280221813315839,0.0354603092047261,-0.0301303275554597,0.378225886411691,-0.143587257286654,0.222998094493725,1.01701617736993,0.0603662733449291,0.526014892794075,0.648051701789525,-0.438204525177852,-0.311563099501274,0.298180173121841,-0.0218836044656445,0.441175053096456,0.357384655770861,0.435621422839663,-0.16604186858177,-0.522763862872567,-0.194570374379559,0.0769622217937904,-0.0322438555767891)
bg188R5 <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df188R5 <- data.frame(exp188R5, bg188R5)
OneWayAnova(data = df188R5, dv = "exp188R5", iv = "bg188R5", conf = 0.95)

#tRF-18-8R6Q46D2
exp188RQ <- c(-0.378059612182765,0.264931809898272,-0.354682819969794,0.0835162816824349,-0.0958618063565827,0.480156146928437,-0.120535855122579,0.130116452826002,0.914910179382668,0.228004663318703,0.489067903223118,0.665286392364235,-0.272572581495748,-0.244635355795665,0.293877348366719,0.0883439209152021,0.363672509558471,0.274990826729237,0.604817328506467,-0.269639251616232,-0.721468304742993,-0.310401771391682,-0.0935979757344345,-0.161756644016547)
bg188RQ <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df188RQ <- data.frame(exp188RQ, bg188RQ)
OneWayAnova(data = df188RQ, dv = "exp188RQ", iv = "bg188RQ", conf = 0.95)

#tRF-18-HR0VX6D2
exp18HRO <- c(-0.272217758348292,-0.0673854547866244,-0.678948064284787,0.206616386507409,0.17236680664456,0.639568084267743,-0.15962446358304,-0.153825115173323,0.615083222294041,0.199180066913279,0.671692323098226,0.69349840820894,-0.360300283209707,-0.581325181426358,-0.212060374875575,-0.0982994262884578,0.50894175806271,0.167206330859478,0.523421789457846,-0.0758289869620246,-0.53218427481942,-0.145221393114039,-0.429057485085791,-0.312076236808442)
bg18HRO <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df18HRO <- data.frame(exp18HRO, bg18HRO)
OneWayAnova(data = df18HRO, dv = "exp18HRO", iv = "bg18HRO", conf = 0.95)

#tRF-228EKSP1852
exp228EK <- c(-0.330512521411821,0.300503426729549,-0.464672096718653,0.088174563325108,0.038471799544359,0.368034828531478,-0.180312602338003,0.249531582439875,1.02183653563478,0.153134163039312,0.608463799502461,0.59196182998028,-0.264159600056402,-0.19350460167292,0.491583409003463,0.114689869984444,0.692924459538544,0.194232491131313,0.436459320169583,-0.0455369608291553,-0.46306229053072,-0.35745197189804,-0.0354461571862556,-0.0428445248339723)
bg228EK <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df228EK <- data.frame(exp228EK, bg228EK)
OneWayAnova(data = df228EK, dv = "exp228EK", iv = "bg228EK", conf = 0.95)

#tRF-22-WE8SPOX52
exp22SPOX <- c(-0.694716029100691,1.57900786185515,-0.506144160762689,-0.314752888593528,-0.335809185591908,0.272414402193678,-0.63748127584164,-0.329147596803988,0.827730174625477,-0.194099734539058,0.106704439891628,0.44486731882728,-0.5602513160078,-0.466569048346557,-0.170291190215541,-0.152488938734459,0.121685310450106,-0.201531173803822,0.46272908964001,-0.666905554745155,-0.943103264282223,-0.728353448543775,-0.393439870592088,-0.622109034470475)
bg22SPOX <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df22SPOX <- data.frame(exp22SPOX, bg22SPOX)
OneWayAnova(data = df22SPOX, dv = "exp22SPOX", iv = "bg22SPOX", conf = 0.95)

#tRF-22-WEKSPM852
exp22WEK <- c(-0.38817669548177,0.0850338426465633,-0.504316876937965,0.222243955769198,0.00136387095074593,0.583851903053229,-0.123107592817252,0.0329933956927634,0.854668169222697,0.262496441206796,0.614318021202846,0.69711800824773,-0.215172383306815,-0.33936224685187,0.305922147380512,0.0687346186439957,0.57492823663973,0.407454044832299,0.586172950202997,0.047563082850198,-0.413867079171398,-0.138272869859456,-0.112828724333403)
bg22WEK <- c("nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","nonstimulated","LPS","LPS","LPS","LPS","LPS","LPS","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","LPS + nicotine","nicotine","nicotine","nicotine","nicotine","nicotine")
df22WEK <- data.frame(exp22WEK, bg22WEK)
OneWayAnova(data = df22WEK, dv = "exp22WEK", iv = "bg22WEK", conf = 0.95)
