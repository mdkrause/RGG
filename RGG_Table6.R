# MET data
dat <- readRDS('./dat/dat.rds')$MET

# crossing nursery
parents <- readRDS('./dat/dat.rds')$parents

# Additive relationship matrix computed with rrBLUP::A.mat()
G <- readRDS('./dat/dat.rds')$G_matrix 
Ginv <- solve(G)

# Loading needed packages
library(dplyr)
library(asreml)
asreml.options(maxit=60,pworkspace='2gb',workspace='1gb')
library(funtimes)

# Models presented in Table 6
EB <- asreml(fixed = eBLUE ~ Y + mappingL, 
             random = ~ G + L + G:L + G:Y + L:Y + G:L:Y,
             family = asr_gaussian(dispersion = 1), 
             weights = weight,
             data = dat)

E0 <- EB # Same as EB, differs in the estimation of RGG

E0V <- asreml(fixed = eBLUE ~ Y + mappingL, 
              random = ~ G + L + G:fa(L,1) + G:idh(Y) + L:Y + G:L:Y,
              family = asr_gaussian(dispersion = 1), 
              weights = weight,
              data = dat)

E0G <- asreml(fixed = eBLUE ~ Y + mappingL, 
              random = ~ vm(G,Ginv) + L + G:L + G:Y + L:Y + G:L:Y,
              family = asr_gaussian(dispersion = 1), 
              weights = weight,
              data = dat)

E0GV <- asreml(fixed = eBLUE ~ Y + mappingL, 
              random = ~ vm(G,Ginv) + L + G:fa(L,1) + G:idh(Y) + L:Y + G:L:Y,
              family = asr_gaussian(dispersion = 1), 
              weights = weight,
              data = dat)

E1 <- asreml(fixed = eBLUE ~ G + Y + mappingL, 
             random = ~ L + G:L + G:Y + L:Y + G:L:Y,
             family = asr_gaussian(dispersion = 1), 
             weights = weight,
             data = dat)

E1V <- asreml(fixed = eBLUE ~ G + Y + mappingL, 
             random = ~ L + G:fa(L,1) + G:idh(Y) + L:Y + G:L:Y,
             family = asr_gaussian(dispersion = 1), 
             weights = weight,
             data = dat)

E2 <- asreml(fixed = eBLUE ~ first_Y_trial:at(check) + Y_num + mappingL, 
             random = ~ L + G:L + G:Y + L:Y + G:L:Y,
             family = asr_gaussian(dispersion = 1), 
             weights = weight,
             data = dat)

E2V <- asreml(fixed = eBLUE ~ first_Y_trial:at(check) + Y_num + mappingL, 
              random = ~ L + G:fa(L,1) + G:idh(Y) + L:Y + G:L:Y,
              family = asr_gaussian(dispersion = 1), 
              weights = weight,
              data = dat)

E3 <- asreml(fixed = eBLUE ~ Y_num + mappingL, 
             random = ~ first_Y_trial:at(check) + L + G:L + G:Y + L:Y + G:L:Y,
             family = asr_gaussian(dispersion = 1), 
             weights = weight,
             data = dat)

E3V <- asreml(fixed = eBLUE ~ Y_num + mappingL, 
              random = ~ first_Y_trial:at(check) + L + G:fa(L,1) + G:idh(Y) + L:Y + G:L:Y,
              family = asr_gaussian(dispersion = 1), 
              weights = weight,
              data = dat)

E4 <- E0 # Same as E0, differs in the estimation of RGG

E4V <- E0V # Same as E0V, differs in the estimation of RGG

E5 <- E1 # Same as E1, differs in the estimation of RGG

E5V <- E1V # Same as E1, differs in the estimation of RGG

E6 <- asreml(fixed = eBLUE ~ check + Y_num + mappingL, 
             random = ~ G + L + Y + G:L + check:L + check:Y,
             data = dat|>subset(trial=='URT'))

# There is a typo in Table 6
# \Sigma_gy was not included in this model
E6V <- asreml(fixed = eBLUE ~ check + Y_num + mappingL, 
             random = ~ G + L + Y + G:fa(L,1) + check:L + check:Y,
             data = dat|>subset(trial=='URT'))

## Extra step:
# Models E7, E7G, and E8, require empirical BLUP values of environments
dat$E <- as.factor(paste0(dat$L, dat$Y))
model <- asreml(eBLUE ~ G,
                random = ~ E + G:E,
                family = asr_gaussian(dispersion = 1),
                weights = weight,
                data = dat|>subset(check=='yes'))
envBlup <- as.matrix(model$coefficients$random)
envBlup <- data.frame(E = levels(dat$E),
                      envBlup = envBlup[match(levels(dat$E),
                                              gsub("E_","",rownames(envBlup)))])
dat <- left_join(dat, envBlup) ; rm(envBlup, model) ; gc()
##

E7 <- asreml(eBLUE ~ envBlup,
             random = ~ G + G:envBlup,
             family = asr_gaussian(dispersion = 1), 
             weights = weight,
             data = dat)

E7G <- asreml(eBLUE ~ envBlup,
              random = ~ vm(G,Ginv) + G:envBlup,
              family = asr_gaussian(dispersion = 1), 
              weights = weight,
              data = dat)

E8 <- asreml(eBLUE ~ envBlup,
             random = ~ first_Y_trial:at(check) + G:envBlup,
             family = asr_gaussian(dispersion = 1), 
             weights = weight,
             data = dat)

## Next step is to estimate RGG

## Extra pieces of code needed
### First year of testing. In this case, the year of PYT
first_test <- dat%>%
  group_by(G)%>%
  filter(check=='no')%>%
  summarise(Y=min(as.numeric(as.character(Y))))%>%
  arrange(Y)
### Second (or last) year of testing. In this case, the year of URT
second_test <- dat%>%
  group_by(G)%>%
  filter(trial=='URT'&check=='no')%>%
  summarise(Y=min(as.numeric(as.character(Y))))%>%
  arrange(Y)
#### Wrapper to get empirical BLUPs or GEBVs
returnBLUP <- function(model, ID = NULL) {
  if(is.null(ID)){
    cat('Provide G names')}
  else{
    tmp <- as.matrix(model$coefficients$random)
    if (length(grep('vm', rownames(tmp))) != 0) {
      rownames(tmp) <-
        gsub("vm\\(G, Ginv)_", "", rownames(tmp), ignore.case = T)
    }
    tmp <- data.frame(G = ID,
                      eBLUP = tmp[match(ID, gsub("G_", "", rownames(tmp)))])
    colnames(tmp)[2] <- paste(substitute(model))
    return(tmp)
  }
}
###

## Indirect measures of RGG

### Getting predicted/estimated values
EB_hat <- returnBLUP(EB, first_test$G)
E0_hat <- returnBLUP(E0, first_test$G)
E0V_hat <- returnBLUP(E0V, first_test$G)
E0G_hat <- returnBLUP(E0G, first_test$G)
E0GV_hat <- returnBLUP(E0GV, first_test$G)
E1_hat <- data.frame(predict.asreml(E1,"G")$pvals)[,c(1:2)] ; colnames(E1_hat)[2] <- 'E1'
E1V_hat <- data.frame(predict.asreml(E1V,"G")$pvals)[,c(1:2)] ; colnames(E1V_hat)[2] <- 'E1V'
E7_hat <- returnBLUP(E7, first_test$G)
E7G_hat <- returnBLUP(E7G, first_test$G)

first_test <- left_join(first_test, EB_hat, by = 'G')
first_test <- left_join(first_test, E0_hat, by = 'G')
first_test <- left_join(first_test, E0V_hat, by = 'G')
first_test <- left_join(first_test, E0G_hat, by = 'G')
first_test <- left_join(first_test, E0GV_hat, by = 'G')
first_test <- left_join(first_test, E1_hat, by = 'G')
first_test <- left_join(first_test, E1V_hat, by = 'G')
first_test <- left_join(first_test, E7_hat, by = 'G')
first_test <- left_join(first_test, E7G_hat, by = 'G')

### Models E4, E4V, E5, and E5V, consider only lines in URT 
E4_hat <- returnBLUP(E4, second_test$G)
E4V_hat <- returnBLUP(E4V, second_test$G)
E5_hat <- E1_hat|>subset(G%in%second_test$G) ; colnames(E5_hat)[2] <- 'E5'
E5V_hat <- E1V_hat|>subset(G%in%second_test$G) ; colnames(E5V_hat)[2] <- 'E5V'

second_test <- left_join(second_test, E4_hat, by = 'G')
second_test <- left_join(second_test, E4V_hat, by = 'G')
second_test <- left_join(second_test, E5_hat, by = 'G')
second_test <- left_join(second_test, E5V_hat, by = 'G')

### Getting true simulated genetic value to compute expected bias from MET
trueGV <- dat%>%group_by(G)%>%
  filter(check=='no')%>%
  summarise(sim_gv=mean(sim_gv))
first_test <- left_join(first_test, trueGV, by = 'G')

### Linear regressions 

#### True trends
modelMET <- lm(sim_gv~Y, first_test) # from MET, explained in "Bias and linearity" 
modelP <- lm(sim_gv~year_crossing_nursery, parents) # from parents, true RGG. Equation 3

### Indirect measures
tmp <- c('G','E0','E0V','E0G','E0GV','sim_gv')

(res1 <- apply(first_test[,-match(tmp,colnames(first_test))],
              2, function(x,Y=first_test$Y) lm(x~Y)$coefficients[2])[-1])

(res2 <- apply(second_test[,-match('G',colnames(second_test))],
               2, function(x,Y=second_test$Y) lm(x~Y)$coefficients[2])[-1])

tmp <- c('Y','E0','E0V','E0G','E0GV')
dat_kappa <- first_test[,match(tmp,colnames(first_test))]
dat_kappa <- dat_kappa%>%group_by(Y)%>%summarise(across(everything(),list(mean))) # Equation 5
colnames(dat_kappa) <- gsub('_1','',colnames(dat_kappa))
(res3 <- apply(dat_kappa, 
               2, function(x,Y=dat_kappa$Y) lm(cumsum(x)~Y)$coefficients[2])[-1])

## Direct measures of RGG
res4 <- c(E2$coefficients$fixed['first_Y_trial:at(check, no)',1],
          E2V$coefficients$fixed['first_Y_trial:at(check, no)',1],
          coef(E3)$random['first_Y_trial:at(check, no)',1],
          coef(E3V)$random['first_Y_trial:at(check, no)',1],
          E6$coefficients$fixed['Y_num',1],
          E6V$coefficients$fixed['Y_num',1],
          coef(E8)$random['first_Y_trial:at(check, no)',1])
names(res4) <- c('E2','E2V','E3','E3V','E6','E6V','E8')

## Linear regression of empirical BLUE values on calendar year
Pheno <- lm(eBLUE ~ Y_num, dat) 

## Summarizing results 
RGG <- data.frame(model = c(names(res1),names(res2),names(res3),names(res4),
                            'Pheno','trueMET','trueRGG'),
                  RGG_hat = c(res1,res2,res3,res4,
                              Pheno$coefficients[2],
                              modelMET$coefficients[2],
                              modelP$coefficients[2]))
rownames(RGG) <- RGG$model
RGG <- 
  RGG[c('EB','E0','E0V','E0G','E0GV',
        'E1','E1V',
        'E2','E2V',
        'E3','E3V',
        'E4','E4V',
        'E5','E5V',
        'E6','E6V',
        'E7','E7G',
        'E8',
        'Pheno',
        'trueMET','trueRGG'),]

RGG$bias <- ((RGG$RGG_hat-RGG['trueRGG',2])/RGG['trueRGG',2])*100 # in percentage
RGG$expected_bias <- ((RGG$RGG_hat-RGG['trueMET',2])/RGG['trueMET',2])*100 # in percentage
RGG[,c(2:4)] <- round(RGG[,c(2:4)],4)
### Keep in mind that in stochastic simulations results are reported as the average of all simulation runs
RGG

# Linearity

## Indirect models
## Sieve-bootstrap version of the Student’s t-test (funtimes). p-values directly reported 
apply(first_test[,-match(c('G','Y'), colnames(first_test))],
               2, function(x) notrend_test(x)$p.value)

apply(second_test[,-match(c('G','Y'), colnames(first_test))],
      2, function(x) notrend_test(x)$p.value)

apply(dat_kappa[,-match('Y', colnames(dat_kappa))],
      2, function(x) notrend_test(x)$p.value)

## Direct models
linearity <- data.frame(model = c('E2','E2V','E3','E3V','E6','E6V','E8'),
                        z_Score = 
                          c(summary(E2,coef=TRUE)$coef.fixed['first_Y_trial:at(check, no)',3],
                            summary(E2V,coef=TRUE)$coef.fixed['first_Y_trial:at(check, no)',3],
                            summary(E3,coef=TRUE)$coef.random['first_Y_trial:at(check, no)',3],
                            summary(E3V,coef=TRUE)$coef.random['first_Y_trial:at(check, no)',3],
                            summary(E6,coef=TRUE)$coef.fixed['Y_num',3],
                            summary(E6V,coef=TRUE)$coef.fixed['Y_num',3],
                            summary(E8,coef=TRUE)$coef.random['first_Y_trial:at(check, no)',3]))

linearity$p_value <- 2*(1-pnorm(linearity$z_Score))
linearity$log_p_value <- -log10(linearity$p_value) # as reported in Fig 7, page 15
number_test <- 1 # Bonferroni correction (S = 225 in the publication)
linearity$significant <- ifelse(linearity$p_value<0.05/number_test,TRUE,FALSE)
linearity
