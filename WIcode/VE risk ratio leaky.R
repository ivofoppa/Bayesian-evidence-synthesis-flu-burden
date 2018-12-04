## Discrete stochastic transmission model with unimdal vaccination uptake, two viruses, all-or-none protection
# bfolder <- 'C:/Users/ifoppa/Documents/GitHub/' ## Define root path (where repository)
bfolder <- 'C:/Users/vor1/Documents/GitHub/' ## Define root path (where repository)
r0 <- 1.8

delta <- 1/4 ## infectious period

seas <- 500 ## Number of days in epidemic
### Generating the vaccination rate over time
vacc <- 0.47 ## cumulative vaccination coverage
phi <- .6

ccratio <- 3 ## control-case ratio
ncrit <- 5 ## Min. required number of cases/controls per stratum
###################################################################################################
##  Seed infections ###############################################################################
###################################################################################################
Ivinit <- 6 ## Number of seed infections
Invinit <- 10 ## Number of seed infections
###################################################################################################
## Setting initial values etc.
fpath <- paste0(bfolder,'Waning-Immunity-artefact/WIcode')
setwd(fpath)

library(binhf)
library(survival)
Ntot <- 2000000

pSI <- function(r0,Iv,Inv){ ## infection probability
  beta <- r0 * delta
  if((beta * (Iv + I2))==0){
    return (c(0,0,1))
  } else {
    lambda <- (beta * (Iv + Inv)) / Ntot
    pinfv <- as.numeric(1-exp(-lambda * phi))
    pinfnv <- as.numeric(1-exp(-lambda))
    return (c(pinfv,pinfnv,1 - pinfv,1 - pinfnv))
  }
}

###################################################################################################
prem <- 1 - exp(-delta) ## Daily removal probability
###################################################################################################
Snvinit <- rbinom(1,Ntot,1 - vacc)
Svinit <- Ntot - Snvinit

Rnvinit <- 0
Rvinit <- 0

inits <- list(Svinit=Svinit,Snvinit=Snvinit,
              Ivinit=Ivinit,Invinit=Invinit,
              Rvinit=Rvinit,Rnvinit=Rnvinit)
###################################################################################################
###################################################################################################
###################################################################################################
Snv <- as.integer(inits$Snvinit)
Sv <- as.integer(inits$Svinit)

Inv <- as.integer(inits$Invinit)
Iv <- as.integer(inits$Ivinit)

Rnv <- as.integer(inits$Rnvinit)
Rv <- as.integer(inits$Rvinit)

time <- 1

## Infections in unvaccinated
pinfnv <-  S1nv/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v)) ### probability that infections with virus 1 are in type 1 individuals
pinf21nv <-  S2nv/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v))  ### probability that infections with virus 1 are in type 2 individuals
pinf01nv <-  S0nv/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v))  ### probability that infections with virus 1 are in type 0 individuals
pinf121nv <-  S12nv/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v))  ### probability that infections with virus 1 are in type 12 individuals
## virus 2
pinf12nv <-  S1nv/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v)) ### probability that infections with virus 1 are in type 1 individuals
pinf22nv <-  S2nv/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v))  ### probability that infections with virus 1 are in type 2 individuals
pinf02nv <-  S0nv/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v))  ### probability that infections with virus 1 are in type 0 individuals
pinf122nv <-  S12nv/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v))  ### probability that infections with virus 1 are in type 12 individuals
## Infections in vaccinated
pinf11v <-  sum(S1v)/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v))
pinf121v <-  1 - pinf11nv - pinf21nv - pinf01nv - pinf121nv - pinf11v  ### probability that infections with virus 1 are in type 1 individuals
## virus 2
pinf22v <-  sum(S2v)/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v))
pinf122v <-  1 - pinf12nv - pinf22nv - pinf02nv - pinf122nv - pinf22v  ### probability that infections with virus 1 are in type 1 individuals

## infection numbers by virus and type
## virus 1
inf1 <- rmultinom(1,inf1num,c(pinf11nv,pinf21nv,pinf01nv,pinf121nv,pinf11v,pinf121v))
inf11nv <- inf1[1]
inf21nv <- inf1[2]
inf01nv <- inf1[3]
inf121nv <- inf1[4]
inf11vnum <- inf1[5]
inf11v <- c(rmultinom(1,inf11vnum,rep(1/prevdur,prevdur)),rep(0,seas))
inf121vnum <- inf1[6]
inf121v <- c(rmultinom(1,inf121vnum,rep(1/prevdur,prevdur)),rep(0,seas))

## virus 2
inf2 <- rmultinom(1,inf2num,c(pinf12nv,pinf22nv,pinf02nv,pinf122nv,pinf22v,pinf122v))
inf12nv <- inf2[1]
inf22nv <- inf2[2]
inf02nv <- inf2[3]
inf122nv <- inf2[4]
inf22vnum <- inf2[5]
inf22v <- c(rmultinom(1,inf22vnum,rep(1/prevdur,prevdur)),rep(0,seas))
inf122vnum <- inf2[6]
inf122v <- c(rmultinom(1,inf122vnum,rep(1/prevdur,prevdur)),rep(0,seas))

#########################################################################################
#########################################################################################
## updating states: Infections with virus 1
I1nv <- I1nv + inf11nv + inf21nv + inf01nv + inf121nv
I1v <- I1v + inf11v + inf121v
## updating states: Infections with virus 2
I2nv <- I2nv + inf12nv + inf22nv + inf02nv + inf122nv
I2v <- I2v + inf22v + inf122v

## Updating # susceptibles
## unvaccinated
S1nv <- S1nv - inf11nv - inf12nv
S2nv <- S2nv - inf21nv - inf22nv
S0nv <- S0nv - inf01nv - inf02nv
S12nv <- S12nv - inf121nv - inf122nv
## vaccinated
S1v <- S1v - inf11v
S2v <- S2v - inf22v
S12v <- S12v - inf121v - inf122v

I1 <- sum(c(I1nv,I1v))
I2 <- sum(c(I2nv,I2v))

I1ls <- I2ls <- kls <- NULL

## Running model
source('model run.R') ## run initialization code

## ... and setting-up data
source('data set prep.R') ## run initialization code
######################################################################################################
######################################################################################################
plot(I1ls, ylim = c(0,max(c(I1ls,I2ls))),type = 'l',col = 'blue')
lines(I2ls,col = 'red')
######################################################################################################
######################################################################################################
# cond_logist <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset,method = 'approximate')
# summary(cond_logist)
# 
cond_logist1 <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset1,method = 'approximate')
summary(cond_logist1)

cond_logist2 <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset12,method = 'approximate')
summary(cond_logist2)

cond_logist3 <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset3,method = 'approximate')
summary(cond_logist3) ## sincevacc not a factor

logist1 <- glm(case ~ sincevacc,weights = count, data = dataset,family = binomial(link = 'logit'))
summary(logist1)

logist2 <- glm(cbind(cases,controls) ~ vacc + late + late*vacc, data = dataset4,family = binomial(link = 'logit'))
summary(logist2)

logist2 <- glm(cbind(cases,controls) ~ vacc + late + late*vacc, data = dataset3,family = binomial(link = 'logit'))
summary(logist2)
######################################################################################################
######################################################################################################

plot(VEls[1:seas2])
lines(trueVEls[1:seas2],col = 'red')

filepath <- paste0('C:/Users/VOR1/Documents/GitHub/Waning-Immunity-artefact/WIwriteup/WIplots/simul_adj_',prevdur,'_',vdur,'.RData')
save(dataset,dataset2,studydata,studydata2,file = filepath)

