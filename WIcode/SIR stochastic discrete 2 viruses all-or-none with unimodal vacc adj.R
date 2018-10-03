## Discrete stochastic transmission model with unimdal vaccination uptake, two viruses, all-or-none protection
library(binhf)
library(survival)
Ntot <- 2000000

r10 <- 1.5
r20 <- 1.6
delta <- 1/4 ## infectious period

prem <- 1 - exp(-delta) ## Daily removal probability

pSE <- function(r0,I){ ## infection probability
  beta <- r0 * delta
  lambda <- beta * I / Ntot
  return (as.numeric(1-exp(-lambda)))
}

seas <- 500 ## Number of days in epidemic
### Generating the vaccination rate over time
vacc1 <- 0.47 ## cumulative vaccination coverage

vs1 <- .3 ## proportion of subject who, if vacc. remain susc. to virus 1
vs2 <- .6 ## proportion of subject who, if vacc. remain to virus 2
vs0 <- 1 - vs1 - vs2

# ## scaling the second virus to about virus 1 transmissability
# fac1 <- vacc1 *vs1 + (1 - vacc1)
# fac2 <- vacc1 *vs2 + (1 - vacc1)
# r20 <- r10 * fac2/fac1

prevdur <- 100 ## vaccination before transmission
vdur <- prevdur + 200 ### duration of vaccination

vdurls <- 1:vdur

raw_vls <- sapply(vdurls,function(x) dnorm(x,vdur/2,vdur/5)) ## Normally distributed vacc uptake
s1 <- -log(vacc1)
v0 <- raw_vls/sum(raw_vls) * s1 ## scaled right
v <- c(v0,rep(0,seas + prevdur - vdur))
vaccdelim <- c(seq(1,vdur,round(vdur/5)),seas + prevdur)
nvaccat <- length(vaccdelim) - 1 ## number of time-since-vacc categories
###################################################################################################
###################################################################################################
S1nvinit <- rbinom(1,Ntot,vs1) ## only susceptible to virus 1 if vacc.
S2nvinit <- rbinom(1,Ntot - S1nvinit,vs2/(1 - vs1)) ## only susceptible to virus 2 if vacc.
S0nvinit <- Ntot - S1nvinit - S2nvinit ## not susceptible to either virus if vacc.

S1vinit <- S2vinit <- S0vinit <-  rep(0,seas + prevdur) ## Susceptible to only virus 1, 2 neither if vacc.

I1nvinit <- I2nvinit <- 0
I1vinit <- I2vinit <-  rep(0,seas + prevdur)

Rnvinit <- 0
Rvinit <- rep(0,seas + prevdur)

inits <- list(S1vinit=S1vinit,S2vinit=S2vinit,S0vinit=S0vinit,
              S1nvinit=S1nvinit,S2nvinit=S2nvinit,S0nvinit=S0nvinit,
              I1vinit=I1vinit,I2vinit=I2vinit,
              I1nvinit=I1nvinit,I2nvinit=I2nvinit,
              Rvinit=Rvinit,Rnvinit=Rnvinit)

# sim <- 1
# ccratio <- 3

S1nv <- as.integer(inits$S1nvinit)
S2nv <- as.integer(inits$S2nvinit)
S0nv <- as.integer(inits$S0nvinit)

S1v <- as.integer(inits$S1vinit)
S2v <- as.integer(inits$S2vinit)
S0v <- as.integer(inits$S0vinit)

I1nv <- as.integer(inits$I1nvinit)
I2nv <- as.integer(inits$I2nvinit)

I1v <- as.integer(inits$I1vinit)
I2v <- as.integer(inits$I2vinit)

Rnv <- as.integer(inits$Rnvinit)

Rv <- as.integer(inits$Rvinit)

time <- 1

while (time <= prevdur){
  ###################################################################################################
  ###################################################################################################
  ### Vaccination
  pvacc <- exp(-v[time]) ## vaccination uptake at this time

  S1v <- shift(S1v,1)
  s1vacc <- rbinom(1,S1nv,1 - pvacc)
  S1v[1] <- s1vacc; S1nv <- S1nv - s1vacc
  
  S2v <- shift(S2v,1)
  s2vacc <- rbinom(1,S2nv,1 - pvacc)
  S2v[1] <- s2vacc; S2nv <- S2nv - s2vacc
  
  S0v <- shift(S0v,1)
  s0vacc <- rbinom(1,S0nv,1 - pvacc)
  S0v[1] <- s0vacc; S0nv <- S0nv - s0vacc
  
  time <- time + 1
}

studydata <- array(0,dim = c(seas,1 + 2 * nvaccat + 2))
studydata2 <- array(0,dim = c(seas,1 + 2 * nvaccat + 2))

inf1num <- 90 ## Number of seed infections
inf2num <- 10 ## Number of seed infections
## Infections in unvaccinated
pinf11nv <-  S1nv/(S1nv + S2nv + S0nv + sum(S1v))  ### probability that infections with virus 1 are in type 1 individuals
pinf21nv <-  S2nv/(S1nv + S2nv + S0nv + sum(S1v))  ### probability that infections with virus 1 are in type 2 individuals
pinf01nv <-  S0nv/(S1nv + S2nv + S0nv + sum(S1v))  ### probability that infections with virus 1 are in type 0 individuals
## virus 2
pinf12nv <-  S1nv/(S1nv + S2nv + S0nv + sum(S2v))  ### probability that infections with virus 1 are in type 1 individuals
pinf22nv <-  S2nv/(S1nv + S2nv + S0nv + sum(S2v))  ### probability that infections with virus 1 are in type 2 individuals
pinf02nv <-  S0nv/(S1nv + S2nv + S0nv + sum(S2v))  ### probability that infections with virus 1 are in type 0 individuals
## Infections in vaccinated
pinf11v <-  1 - pinf11nv - pinf21nv - pinf01nv  ### probability that infections with virus 1 are in type 1 individuals
## virus 2
pinf22v <-  1 - pinf11nv - pinf21nv - pinf01nv  ### probability that infections with virus 1 are in type 1 individuals

## infection numbers by virus and type
## virus 1
inf11nv <- rbinom(1,inf1num,pinf11nv)
inf21nv <- rbinom(1,inf1num - inf11nv, pinf21nv/(pinf21nv + pinf01nv + pinf11v))
inf01nv <- rbinom(1,inf1num - inf11nv - inf21nv, pinf01nv/(pinf01nv + pinf11v))
inf11vnum <- inf1num - inf11nv - inf21nv - inf01nv
inf11v <- c(rmultinom(1,inf11vnum,rep(1/prevdur,prevdur)),rep(0,seas)); S1v <- S1v - inf11v

## virus 2
inf12nv <- rbinom(1,inf2num,pinf12nv)
inf22nv <- rbinom(1,inf2num - inf12nv, pinf22nv/(pinf22nv + pinf02nv + pinf22v))
inf02nv <- rbinom(1,inf2num - inf12nv - inf22nv, pinf02nv/(pinf02nv + pinf22v))
inf22vnum <- inf2num - inf12nv - inf22nv - inf02nv
inf22v <- c(rmultinom(1,inf22vnum,rep(1/prevdur,prevdur)),rep(0,seas)); S2v <- S2v - inf22v
#########################################################################################
#########################################################################################
## updating states: Infections with virus 1
I1nv <- I1nv + inf11nv + inf21nv + inf01nv
I1v <- I1v + inf11v
## updating states: Infections with virus 2
I2nv <- I2nv + inf12nv + inf22nv + inf02nv
I2v <- I2v + inf22v
## updating states: Infections with either virus
## unvaccinated
S1nv <- S1nv - inf11nv - inf12nv
S2nv <- S2nv - inf21nv - inf22nv
S0nv <- S0nv - inf01nv - inf02nv
## vaccinated
S1v <- S1v - inf11v
S2v <- S2v - inf22v

I1 <- sum(c(I1nv,I1v))
I2 <- sum(c(I2nv,I2v))

I1ls <- I2ls <- NULL

while (time <= seas + prevdur){
  ## virus 1
  p1inf <- pSE(r10,I1)
  inf11nv <- rbinom(1, S1nv, p1inf)
  inf21nv <- rbinom(1, S2nv, p1inf)
  inf01nv <- rbinom(1, S0nv, p1inf)
  
  inf11v <- sapply(S1v, function(sv) rbinom(1, sv, p1inf))
  inf11v2 <- sapply(seq_along(inf11v), function(k) ifelse(k <= (time - prevdur), 0, inf11v[k]))
  
  rem1nv <- rbinom(1,I1nv,prem)
  rem1v <- sapply(I1v, function(iv) rbinom(1,iv, prem))
  ## for plotting purposes ...
  I1ls <- c(I1ls,sum(inf11v) + inf11nv + inf21nv + inf01nv)  
  ## virus 2
  p2inf <- pSE(r20,I2)
  inf12nv <- rbinom(1, S1nv, p2inf)
  inf22nv <- rbinom(1, S2nv, p2inf)
  inf02nv <- rbinom(1, S0nv, p2inf)
  
  inf22v <- sapply(S1v, function(sv) rbinom(1, sv, p2inf))
  inf22v2 <- sapply(seq_along(inf22v), function(k) ifelse(k <= (time - prevdur), 0, inf22v[k]))
  
  rem2nv <- rbinom(1,I2nv,prem)
  rem2v <- sapply(I2v, function(iv) rbinom(1,iv, prem))
  ## for plotting purposes ...
  I2ls <- c(I2ls,sum(inf22v) + inf12nv + inf22nv + inf02nv)  
  #########################################################################################
  #########################################################################################
  ## updating states: Infections with virus 1
  I1nv <- I1nv + inf11nv + inf21nv + inf01nv - rem1nv
  I1v <- I1v + inf11v - rem1v
  ## updating states: Infections with virus 2
  I2nv <- I2nv + inf12nv + inf22nv + inf02nv - rem2nv
  I2v <- I2v + inf22v - rem2v
  ## updating states: Infections with either virus
  ## unvaccinated
  S1nv <- S1nv - inf11nv - inf12nv
  S2nv <- S2nv - inf21nv - inf22nv
  S0nv <- S0nv - inf01nv - inf02nv
  Rnv <- Rnv + rem1nv + rem2nv
  ## vaccinated
  S1v <- S1v - inf11v
  S2v <- S2v - inf22v
  Rv <- Rv + rem1v + rem2v
  ## Total numbers infectious with the two viruses  
  I1 <- sum(c(I1nv,I1v))
  I2 <- sum(c(I2nv,I2v))
  ###################################################################################################
  ### "Study" is conducted ##########################################################################
  ###################################################################################################
  infnv <- inf11nv + inf21nv + inf01nv + inf12nv + inf22nv + inf02nv
  infv <- inf11v + inf22v
  infv2 <- inf11v2 + inf22v2
  
  casesls <- c(infnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(infv[vaccdelim[x] : vaccdelim[x + 1]])))
  casesls2 <- c(infnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(infv2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  controlsvsum <- S1v + S2v + S0v + Rv
  controlsvsum2 <- sapply(seq_along(controlsvsum), function(k) ifelse(k <= (time - prevdur), 0, controlsvsum[k]))
  
  controlsnv <- S1nv + S2nv + S0nv + Rnv
  
  controlsls <- c(controlsnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum[vaccdelim[x] : vaccdelim[x + 1]])))
  controlsls2 <- c(controlsnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  studydata[time - prevdur,] <- c(time - prevdur,casesls,controlsls)
  studydata2[time - prevdur,] <- c(time - prevdur,casesls2,controlsls2)
  ###################################################################################################
  ### Stop when no more transmission 
  if (I1==0 & I2==0){
    break()
  }
  ###################################################################################################
  ### Vaccination: proceeding time since vacc. and adding new vaccinees
  pvacc <- exp(-v[time]) ## vaccination uptake at this time
  
  S1v <- shift(S1v,1)
  s1vacc <- rbinom(1,S1nv,1 - pvacc)
  S1v[1] <- s1vacc; S1nv <- S1nv - s1vacc
  
  S2v <- shift(S2v,1)
  s2vacc <- rbinom(1,S2nv,1 - pvacc)
  S2v[1] <- s2vacc; S2nv <- S2nv - s2vacc
  
  S0v <- shift(S0v,1)
  s0vacc <- rbinom(1,S0nv,1 - pvacc)
  S0v[1] <- s0vacc; S0nv <- S0nv - s0vacc
  
  I1v <- shift(I1v,1) ### Infectious not getting vaccinated, by assumption
  I2v <- shift(I2v,1) ### 
  # 
  Rv <- shift(Rv,1)

  time <- time + 1
}

delind <- which(rowSums(studydata)==0 | any(is.na(studydata)))

studydata <- studydata[-delind,]
studydata2 <- studydata2[-delind,]

seas2 <- head(delind,1) - 1
### Reorganizing data set for analysis--only vaccinated
dataset <- dataset2 <- NULL

for (t in 1:seas2){
  # totcases <- sum(studydata[t,2:7])
  # totnoncases <- sum(studydata[t,8:13])
  # oddscontrol <- totcases/totnoncases * ccratio
  # pcontrol <- oddscontrol * totnoncases/Ntot 
  # 
  # totcases2 <- sum(studydata2[t,2:7])
  # totnoncases2 <- sum(studydata2[t,8:13])
  # oddscontrol2 <- totcases2/totnoncases2 * ccratio
  # pcontrol2 <- oddscontrol2 * totnoncases2/Ntot 
  
  pcontrol <- pcontrol2 <- 1
  for (k in seq_along(vaccdelim[-1])){
    ncases <- studydata[t,2 + k]
    nnoncases <- studydata[t,8 + k]
    
    ncases2 <- studydata2[t,2 + k]
    nnoncases2 <- studydata2[t,8 + k]
    
    # ncontrols <- rbinom(1,nnoncases,pcontrol)
    # ncontrols2 <- rbinom(1,nnoncases2,pcontrol2)
    ncontrols <- nnoncases
    ncontrols2 <- nnoncases2
    if ((ncases > 5 & ncontrols > 5) & (!is.na(ncases) & !is.na(ncontrols))){
      datarec <- c(t,k,1,ncases)
      dataset <- rbind(dataset,datarec,deparse.level = 0)
      datarec <- c(t,k,0,ncontrols)
      dataset <- rbind(dataset,datarec,deparse.level = 0)
    }
    if ((ncases2 > 5 & ncontrols2 > 5) & (!is.na(ncases2) & !is.na(ncontrols2))){
      datarec2 <- c(t,k,1,ncases2)
      dataset2 <- rbind(dataset2,datarec2,deparse.level = 0)
      datarec2 <- c(t,k,0,ncontrols2)
      dataset2 <- rbind(dataset2,datarec2,deparse.level = 0)
    }
  }
}

colnames(dataset) <- c('time','sincevacc','case','count')
dataset <- data.frame(dataset)
dataset$sincevacc <- factor(dataset$sincevacc)

colnames(dataset2) <- c('time','sincevacc','case','count')
dataset2 <- data.frame(dataset2)
dataset2$sincevacc <- factor(dataset2$sincevacc)

# delind <- which(dataset$case==1 & dataset$count==0)
# dataset <- dataset[-delind,]

cond_logist <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset,method = 'approximate')
summary(cond_logist)

cond_logist2 <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset2,method = 'approximate')
summary(cond_logist2)

<<<<<<< HEAD
# filepath <- paste0('C:/Users/VOR1/Documents/GitHub/Waning-Immunity-artefact/WIwriteup/WIplots/','simul_2_virus_adj_100_300.RData')
# load(filepath)
filepath <- paste0('C:/Users/VOR1/Documents/GitHub/Waning-Immunity-artefact/WIwriteup/WIplots/simul_2_virus_adj_',prevdur,'_',vdur,'.RData')
=======
filepath <- paste0('C:/Users/VOR1/Documents/GitHub/Waning-Immunity-artefact/WIwriteup/WIplots/simul_adj_',prevdur,'_',vdur,'.RData')
>>>>>>> parent of 8315afe... update
save(dataset,dataset2,studydata,studydata2,file = filepath)
######################################################################################################
##  Crude analysis, only adjusting for time ##########################################################
##  Create new data set first               ##########################################################
######################################################################################################
dataset3 <- NULL
VEls <- NULL
ccratio <- 2
for (t in 1:seas2){
  totcases <- sum(studydata[t,2:7])
  totnoncases <- sum(studydata[t,8:13])
  oddscontrol <- totcases/totnoncases * ccratio
  pcontrol <- oddscontrol * totnoncases/Ntot
  
  ncasesnv <- studydata[t,2]
  nnoncasesnv <- studydata[t,8]
  
  ncasesv <- sum(studydata[t,3:7])
  nnoncasesv <- sum(studydata[t,9:13])
  
  ncontrolsv <- rbinom(1,nnoncasesv,pcontrol)
  ncontrolsnv <- rbinom(1,nnoncasesnv,pcontrol)
  # ncontrolsv <- nnoncasesv
  # ncontrolsnv <- nnoncasesnv
  
  if (ncasesv >= 5 & ncasesnv >= 5){
    datarec <- c(t,1,ncasesv,ncontrolsv)
    dataset3 <- rbind(dataset3,datarec,deparse.level = 0)
    datarec <- c(t,0,ncasesnv,ncontrolsnv)
    dataset3 <- rbind(dataset3,datarec,deparse.level = 0)
    
    veest <- 1 - ncasesv * ncontrolsnv / (ncasesnv * ncontrolsv)
    VEls <- c(VEls, veest)
  }
}
colnames(dataset3) <- c('time','vacc','cases','controls')
dataset3 <- data.frame(dataset3)

logist <- glm(cbind(cases,controls) ~ vacc, data = dataset3,family = binomial(link = 'logit'))
summary(logist)
######################################################################################################
######################################################################################################
plot(VEls)

plot(I1ls)
lines(I2ls)
