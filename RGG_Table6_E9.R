## Note: the code below can be optimized in several ways 

# MET data
dat <- readRDS('dat_seed_45.rds')$MET|>subset(trial!='BT3')

# crossing nursery
parents <- readRDS('dat_seed_45.rds')$parents

# Loading needed packages
library(asreml)
asreml.options(maxit=60,pworkspace='2gb',workspace='1gb',trace=FALSE)
library(dplyr)
library(tidyr)

# All years are used to correct the data
years <- levels(dat$Y)
eblup <- c()

# analysis within years to estimate Cullins H2 and compute eBLUPs
for(i in years){
  datTemp <- dat %>% filter(Y == i) %>% droplevels()
  model <- asreml(eBLUE ~ -1 + L,  
                  random = ~G + L:G, 
                  family = asr_gaussian(dispersion = 1), 
                  weights = weight,
                  data = datTemp)
  bT <- as.matrix(model$coefficients$random)                       
  bT <- data.frame(G = levels(datTemp$G), 
                   bT = bT[match(levels(datTemp$G),gsub("G_","",rownames(bT)))])
  bT$intercept <- mean(model$coefficients$fixed[,1])
  bT$Y <- i
  
  vdBLUP.mat <- predict.asreml(model, classify="G", only="G", sed=TRUE)$sed^2  
  vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)])
  vg <- summary(model)$varcomp['G', 'component']
  H2Cullis <- 1 - (vdBLUP.avg/(vg*2))
  bT$H2 <- H2Cullis
  eblup <- rbind(eblup, bT)
  rm(datTemp, bT, model, vdBLUP.mat, vdBLUP.avg, H2Cullis, vg)
  print(i)
}

# Choosing reference year according to H2
eblup$pheno <- eblup$bT+eblup$intercept
eblup$Y <- as.numeric(as.character(eblup$Y))

H2 <- eblup %>% group_by(Y) %>% summarise(H2 = unique(H2))
(H2 <- H2 %>% top_n(5)) # top 5
(ref <- sample(H2$Y,1)) # sample the reference year
H2 <- H2 %>% filter(Y == ref) 

refInit <- ref
#checks <- dat %>% group_by(G) %>% filter(check == 'no' & trial == 'URT') %>% summarise(Y = max(as.numeric(as.character(Y))))
#checks$Y <- as.character(checks$Y)

# algorithm for correcting the "year effect" according to the reference year

# forward
if(ref == min(eblup$Y)){
  for(ref in ref:(max(eblup$Y)-1)){
    common <- eblup %>% filter(Y %in% c(ref, ref+1)) %>% droplevels()
    tab <- table(common$G, common$Y) 
    tab <- apply(tab, 1, sum)
    common <- common %>% filter(G %in% names(which(tab > 1)))
    common <- common %>% select(G, Y, pheno) %>% spread(key = Y, value = pheno)
    correc <- sum(common[,3]-common[,2])/nrow(common)
    eblup$pheno[which(eblup$Y == c(ref+1))] <- eblup$pheno[which(eblup$Y == c(ref+1))]+correc
    #print(ref)
    rm(common, tab, correc)
  }
}

# backward
if(ref == max(eblup$Y)){
  for(ref in ref:(min(eblup$Y)+1)){
    common <- eblup %>% filter(Y %in% c(ref, ref-1)) %>% droplevels()
    tab <- table(common$G, common$Y) 
    tab <- apply(tab, 1, sum)
    common <- common %>% filter(G %in% names(which(tab > 1)))
    common <- common %>% select(G, Y, pheno) %>% spread(key = Y, value = pheno)
    correc <- sum(common[,2]-common[,3])/nrow(common)
    eblup$pheno[which(eblup$Y == c(ref-1))] <- eblup$pheno[which(eblup$Y == c(ref-1))]+correc
    rm(common, tab, correc)
  }
}


# forward-backward
if(between(ref, min(eblup$Y)+1, max(eblup$Y)-1)){
  
  # forward
  for(ref in ref:(max(eblup$Y)-1)){
    common <- eblup %>% filter(Y %in% c(ref, ref+1)) %>% droplevels()
    tab <- table(common$G, common$Y) 
    tab <- apply(tab, 1, sum)
    common <- common %>% filter(G %in% names(which(tab > 1)))
    common <- common %>% select(G, Y, pheno) %>% spread(key = Y, value = pheno)
    correc <- sum(common[,3]-common[,2])/nrow(common)
    eblup$pheno[which(eblup$Y == c(ref+1))] <- eblup$pheno[which(eblup$Y == c(ref+1))]+correc
    rm(common, tab, correc)
  }
  
  ref <- refInit
  
  #backward
  for(ref in ref:(min(eblup$Y)+1)){
    common <- eblup %>% filter(Y %in% c(ref, ref-1)) %>% droplevels()
    tab <- table(common$G, common$Y) 
    tab <- apply(tab, 1, sum)
    common <- common %>% filter(G %in% names(which(tab > 1)))
    common <- common %>% select(G, Y, pheno) %>% spread(key = Y, value = pheno)
    correc <- sum(common[,2]-common[,3])/nrow(common)
    eblup$pheno[which(eblup$Y == c(ref-1))] <- eblup$pheno[which(eblup$Y == c(ref-1))]+correc
    rm(common, tab, correc)
  }
}

# the corrections ends here. Let's compute RGG:

corrected_dat <- dat %>% filter(check == 'no') %>% distinct(G, Y) # only URT
#corrected_dat <- dat %>% filter(check == 'no' & trial == 'URT') %>% distinct(G, Y) # PYT + URT
corrected_dat$Y <- as.numeric(as.character(corrected_dat$Y))
corrected_dat <- left_join(corrected_dat, eblup[,c('G', 'Y', 'pheno')])

# True trend
modelP <- lm(sim_gv~year_crossing_nursery, parents) # from parents, true RGG. Equation 3

# RGG from model E9
E9 <- lm(pheno~Y, corrected_dat) 

# Summarizing results
RGG <- data.frame(trueRGG = modelP$coefficients[2],
                   E9 =  E9$coefficients[2],
                   refYear = refInit,
                   H2 = H2$H2) ; rownames(RGG) <- NULL
RGG
