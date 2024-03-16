 # MET data
dat <- readRDS('dat_seed_45.rds')$MET

# crossing nursery
parents <- readRDS('dat_seed_45.rds')$parents

# Additive relationship matrix computed with rrBLUP::A.mat()
G <- readRDS('dat_seed_45.rds')$G_table7
Ginv <- solve(G)

# Loading needed packages
library(dplyr)
library(asreml)
asreml.options(maxit=60,pworkspace='7gb',workspace='5gb')
library(funtimes)

# adding populations, which in the paper I called "cross"
dat$cross <- as.factor(paste0(dat$P1, dat$P2))

# modeling
E10 <- asreml(fixed = eBLUE ~ Y + mappingL, 
             random = ~ P1 + and(P2) + cross + L + L:Y + 
               P1:L + P2:L + cross:L +
               P1:Y + P2:Y + cross:Y +
               P1:L:Y + P2:L:Y,
             equate.levels=c('P1','P2'),
             residual=~idv(units),
             data = dat)

E10G <- asreml(fixed = eBLUE ~ Y + mappingL, 
              random = ~ vm(P1, Ginv) + and(vm(P2, Ginv)) + cross + L + L:Y + 
                P1:L + P2:L + cross:L +
                P1:Y + P2:Y + cross:Y +
                P1:L:Y + P2:L:Y,
              equate.levels=c('P1','P2'),
              residual=~idv(units),
              data = dat)

E11 <- asreml(fixed = eBLUE ~ Y + mappingL + P1 + and(P2) + cross,
              random = ~ L + L:Y +
                P1:L + P2:L + cross:L +
                P1:Y + P2:Y + cross:Y +
                P1:L:Y + P2:L:Y +
                cross:L:Y,
              equate.levels=c('P1','P2'),
              residual=~idv(units),
              data = dat)

E12 <- asreml(fixed = eBLUE ~ Y + mappingL + first_Y_trial:at(parent),
              random = ~ cross + L + L:Y +
                P1:L + P2:L + cross:L +
                P1:Y + P2:Y + cross:Y +
                P1:L:Y + P2:L:Y +
                cross:L:Y,
              equate.levels=c('P1','P2'),
              residual=~idv(units),
              data = dat)
E12 <- update.asreml(E12)

## Extra step for Models E13 and E14
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

E13 <- asreml(fixed = eBLUE ~ envBlup + P1 + and(P2),
              random = ~ cross + P1:envBlup + P2:envBlup + cross:envBlup,
              equate.levels=c('P1','P2'),
              residual=~idv(units),
              data = dat)

E14 <- asreml(fixed = eBLUE ~ envBlup,
              random = ~ first_Y_trial:at(parent) + cross + cross:envBlup,
              residual=~idv(units),
              data = dat)

## Next step is to estimate RGG

## Extra pieces of code needed
### First year of testing for parents/breeding lines
P1 <- dat %>% filter(check == 'no')%>% group_by(P1) %>% summarise(Y = min(as.numeric(as.character(Y))))
P2 <- dat  %>% filter(check == 'no') %>% group_by(P2) %>% summarise(Y = min(as.numeric(as.character(Y))))
colnames(P1)[1] <- 'G' ; colnames(P2)[1] <- 'G'
first_test_P <- rbind(P1, P2) ; rm(P1,P2)
first_test_P <- first_test_P %>% group_by(G) %>% arrange(Y) %>% distinct() 
#### Wrapper to get empirical BLUPs or GEBVs
returnBLUP <- function(model, ID = NULL) {
  if(is.null(ID)){
    cat('Provide G names')}
  else{
    tmp <- as.matrix(model$coefficients$random)
    if (length(grep('vm', rownames(tmp))) != 0) {
      rownames(tmp) <-
        gsub("vm\\(P1, Ginv)_", "", rownames(tmp), ignore.case = T)
    }
    tmp <- data.frame(G = ID,
                      eBLUP = tmp[match(ID, gsub("P1_", "", rownames(tmp)))]*2)
    colnames(tmp)[2] <- paste(substitute(model))
    return(tmp)
  }
}
###

## Indirect measures of RGG

### Getting predicted/estimated values
E10_hat <- returnBLUP(E10, first_test_P$G)
E10G_hat <- returnBLUP(E10G, first_test_P$G)

E11_hat <- as.matrix(coef(E11)$fixed)                       
E11_hat <- data.frame(G = first_test_P$G, 
                   E11 = E11_hat[match(first_test_P$G,
                                       gsub("P1_","",rownames(E11_hat)))]*2) %>% na.omit()

E13_hat <- as.matrix(coef(E13)$fixed)                       
E13_hat <- data.frame(G = first_test_P$G, 
                      E13 = E13_hat[match(first_test_P$G,
                                          gsub("P1_","",rownames(E13_hat)))]*2) %>% na.omit()

first_test_P <- left_join(first_test_P, E10_hat, by = 'G')
first_test_P <- left_join(first_test_P, E11_hat, by = 'G')
first_test_P <- left_join(first_test_P, E13_hat, by = 'G')

### Getting true simulated genetic valueS to compute expected bias from MEs
# Note the first few parents won't be present in MET due to the structure of the simulation
trueGV <- dat%>%group_by(G)%>%
  filter(parent=='yes')%>%
  summarise(sim_gv=mean(sim_gv))
first_test_P <- left_join(first_test_P, trueGV, by = 'G')

### Linear regressions 

#### True trends
modelMET <- lm(sim_gv~Y, first_test_P) # from MET, explained in "Bias and linearity".
                                       # Given the models in Table 7 estimate RGG only with
                                       # parents, the slopes from modelMET and modelP are 
                                       # very similar, as expected. 
                                       # The numerical differences between these slopes are due to
                                       # the first few parents not being included in MET
modelP <- lm(sim_gv~year_crossing_nursery, parents) # from parents, true RGG. Equation 3

### Indirect measures

tmp <- c('Y','E10','E11','E13')
dat_kappa <- first_test_P[,match(tmp,colnames(first_test_P))]
dat_kappa <- dat_kappa%>%group_by(Y)%>%summarise(across(everything(),list(mean))) # Equation 5
colnames(dat_kappa) <- gsub('_1','',colnames(dat_kappa))
(res1 <- apply(dat_kappa, 
               2, function(x,Y=dat_kappa$Y) lm(cumsum(x)~Y)$coefficients[2])[-1])

## Direct measures of RGG

res2 <- c(E12$coefficients$fixed['first_Y_trial:at(parent, yes)',1],
          coef(E14)$random['first_Y_trial:at(parent, yes)',1])
names(res2) <- c('E12','E14')

## Summarizing results 

RGG <- data.frame(model = c(names(res1),names(res2),
                            'trueMET','trueRGG'),
                  RGG_hat = c(res1,res2,
                              modelMET$coefficients[2],
                              modelP$coefficients[2]))
rownames(RGG) <- RGG$model
RGG <- 
  RGG[c('E10','E11','E12','E13','E14',
        'trueMET','trueRGG'),]

RGG$bias <- ((RGG$RGG_hat-RGG['trueRGG',2])/RGG['trueRGG',2])*100 # in percentage
RGG$expected_bias <- ((RGG$RGG_hat-RGG['trueMET',2])/RGG['trueMET',2])*100 # in percentage
rownames(RGG) <- NULL
RGG[,c(2:4)] <- round(RGG[,c(2:4)],4) ; RGG
### Keep in mind that in stochastic simulations results are reported 
### as the average of all simulation runs. Here I am presenting only one rep.

# For linearity metrics, please check "RGG_Table6"
