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
  if((beta * (Iv + Inv))==0){
    return (c(0,0,1))
  } else {
    lambda <- (beta * (Iv + Inv)) / Ntot
    pinfv <- as.numeric(1-exp(-lambda * phi))
    pinfnv <- as.numeric(1-exp(-lambda))
    return (c(pinfv,pinfnv))
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

## Running model
studydata1 <- studydata2 <- NULL

while (time <= seas){
  pinf <- pSI(r0,Iv,Inv) ## for fully susceptible
  ## unvaccinated
  infv <- rbinom(1, Sv, pinf[1])
  ## vaccinated
  infnv <- rbinom(1, Snv, pinf[2])

  ### updating susceptibles
  Sv <- Sv - infv
  Snv <- Snv - infnv
  ## updating other states: Infections
  Iv <- Iv + infv
  Inv <- Inv + infnv

  remv <- rbinom(1,Iv, prem)
  remnv <- rbinom(1,Inv, prem)
  ## numbers removed from infectous
  Iv <- Iv - remv
  Inv <- Inv - remnv
  ## numbers removed from virus infectious
  Rv <- Rv + remv
  Rnv <- Rnv + remnv
  
  ### Interrupt evaluation if no more transmission
  if (Iv==0 & Inv==0 | is.na(Iv) | is.na(Inv)){
    break()
  }
  ###################################################################################################
  ### "Study" is conducted ##########################################################################
  ###################################################################################################
  ncases <- infv + infnv
  ncontrols <- ccratio*ncases
  ncontrolel <- Ntot - Iv - Inv
  
  sratio <- ncontrols/ncontrolel
  controlsv1 <- rbinom(1,Sv + Rv,sratio)
  controlsnv1 <- rbinom(1,Snv + Rnv,sratio)

  controlsv2 <- rbinom(1,Svinit,sratio)
  controlsnv2 <- rbinom(1,Snvinit,sratio)
  
  datarr1 <- rbind(c(time,1,1,infv),c(time,1,0,infnv),c(time,0,1,controlsv1),c(time,0,0,controlsnv1),deparse.level = 0)
  datarr2 <- rbind(c(time,1,1,infv),c(time,1,0,infnv),c(time,0,1,controlsv2),c(time,0,0,controlsnv2),deparse.level = 0)
  
  studydata1 <- rbind(studydata1,datarr1,deparse.level = 0)
  studydata2 <- rbind(studydata2,datarr2,deparse.level = 0)
  ###################################################################################################
  ###################################################################################################
  ### Stop when no more transmission 
  time <- time + 1
}

studydata1 <- data.frame(studydata1)
studydata2 <- data.frame(studydata2)
colnames(studydata1) <- colnames(studydata2) <- c('day','case','vacc','count')
######################################################################################################
######################################################################################################
######################################################################################################
logist1 <- glm(case ~ vacc,weights = count, data = studydata1,family = binomial())
logist2 <- glm(case ~ vacc,weights = count, data = studydata2,family = binomial())

c1 <- sum(studydata0$count[which(studydata0$case==1 & studydata0$vacc==1)])
c0 <- sum(studydata0$count[which(studydata0$case==1 & studydata0$vacc==0)])

1 - c1/Svinit/(c0/Snvinit)
1-exp(logist1$coefficients[2])
1-exp(logist2$coefficients[2])
