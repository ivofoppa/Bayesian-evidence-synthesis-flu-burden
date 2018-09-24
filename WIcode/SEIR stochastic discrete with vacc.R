## Discrete stochastic transmission model with contact matrix
## Latent period (from Comflu_bas): 30% 1 day; 50% 2 days and 20% 3 days;
## Infectious period: 30% 3, 40% 4, 20% 5, 10% 6 days
library(binhf)
library(survival)
Ntot <- 2000000
## Creating population according to latent (rows) and infectious periods (columns) in each of the age groups
r0 <- 1.8
gamma <- 1/2 ## latent period
delta <- 1/4 ## infectious period
beta <- r0 * delta

phi <- .6

pSE <- function(r0,I,phi,vacc){ ## infection probability
  lambda <- beta * (1-phi)^vacc * I / Ntot
  return (as.numeric(1-exp(-lambda)))
}

seas <- 500 ## Number of days in epidemic
### Generating the vaccination rate over time
vacc <- .5

prevdur <- 60 ## vaccination before transmission
vdur <- prevdur + 150 ### duration of vaccination

vrate <- -log(1 - vacc)/vdur

v <- c(rep(vrate,vdur),rep(0,seas - vdur))
vaccdelim <- c(seq(1,vdur,round(vdur/5)),seas + prevdur)
nvaccat <- length(vaccdelim) - 1 ## number of time-since-vacc categories
###################################################################################################
###################################################################################################
Invinit <- 0
Envinit <- 0
Snvinit <- Ntot
Rnvinit <- 0

Evinit <- rep(0,seas + prevdur)
Ivinit <- rep(0,seas + prevdur)
Svinit <- rep(0,seas + prevdur)
Rvinit <- rep(0,seas + prevdur)

inits <- list(Svinit=Svinit,Snvinit=Snvinit,Evinit=Evinit,Envinit=Envinit,
              Ivinit=Ivinit,Invinit=Invinit,Rvinit=Rvinit,Rnvinit=Rnvinit)

# sim <- 1
ccratio <- 3

Snv <- as.integer(inits$Snvinit)
Env <- as.integer(inits$Envinit)
Inv <- as.integer(inits$Invinit)
Rnv <- as.integer(inits$Rnvinit)

Sv <- as.integer(inits$Svinit)
Ev <- as.integer(inits$Evinit)
Iv <- as.integer(inits$Ivinit)
Rv <- as.integer(inits$Rvinit)

time <- 1
while (time <= prevdur){
  ###################################################################################################
  ###################################################################################################
  ### Vaccination
  Sv <- shift(Sv,1)
  newvacc <- rbinom(1,Snv,1 - exp(-v[time]))
  Sv[1] <- newvacc; Snv <- Snv - newvacc
  
  time <- time + 1
}

studydata <- array(0,dim = c(seas,1 + 2 * nvaccat + 2))

infnum <- 100

pinfnv <-  Snv/(Snv + (1 - phi)*sum(Sv))  ### probability that infections are non-vaccinated
infnv <- rbinom(1,infnum,pinfnv)
Env <- Env + infnv

infv <- infnum - infnv
pSv <- Sv/max(sum(Sv),1)
Ev <- Ev + c(rmultinom(1,infv,pSv))

I <- sum(c(Iv,Inv))

Ils <- NULL

while (time <= seas + prevdur){
  pnvinf <- pSE(r0,I,phi,0)
  pvinf <- pSE(r0,I,phi,1)
  
  newlatnv <- rbinom(1, Snv, pnvinf)
  newlatv <- sapply(Sv, function(sv) rbinom(1, sv, pvinf))
  
  newinfnv <- rbinom(1,Env,1 - exp(-gamma))
  newinfv <- sapply(Ev,function(ev) rbinom(1,ev,1 - exp(-gamma)))
  
  newremnv <- rbinom(1,Inv,1 - exp(-delta))
  newremv <- sapply(Iv, function(iv) rbinom(1,iv, 1 - exp(-delta)))

  Ils <- c(Ils,sum(newinfv) + newinfnv)  
  Snv <- Snv - newlatnv
  Sv <- Sv - newlatv
  
  Env <- Env + newlatnv - newinfnv
  Ev <- Ev + newlatv - newinfv
  
  Inv <- Inv + newinfnv - newremnv
  Iv <- Iv + newinfv - newremv
  
  I <- sum(c(Iv,Inv))
  
  Rnv <- Rnv + newremnv
  Rv <- Rv + newremv
  
  ###################################################################################################
  ### "Study" is conducted ##########################################################################
  ###################################################################################################
  casesls <- c(newinfnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(newinfv[vaccdelim[x] : vaccdelim[x + 1]])))
  
  controlsvsum <- Sv + Ev + Rv
  controlsls <- c(Snv + Env + Rnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum[vaccdelim[x] : vaccdelim[x + 1]])))
  
  studydata[time - prevdur,] <- c(time - prevdur,casesls,controlsls)
  ###################################################################################################
  ### Stop when no more transmission 
  if (I==0 & sum(c(Env + Ev)) == 0){
    break()
  }
  ###################################################################################################
  ### Vaccination: proceeding time since vacc. and adding new vaccinees
  Sv <- shift(Sv,1)
  svacc <- rbinom(1,Snv,1 - exp(-v[time]))
  Sv[1] <- svacc; Snv <- Snv - svacc
  
  Ev <- shift(Ev,1)
  evacc <- rbinom(1,Env,1 - exp(-v[time]))
  Ev[1] <- evacc; Env <- Env - evacc
  
  Iv <- shift(Iv,1) ### Infectious not vaccinated, by assumption
  # ivacc <- rbinom(1,Inv,1 - exp(-v[time]))
  # Iv[1] <- ivacc; Inv <- Inv - ivacc
  # 
  Rv <- shift(Rv,1)
  rvacc <- rbinom(1,Rnv,1 - exp(-v[time]))
  Rv[1] <- rvacc; Rnv <- Rnv - rvacc
  
  time <- time + 1
}

delind <- which(rowSums(studydata)==0 | any(is.na(studydata)))

studydata <- studydata[-delind,]

seas2 <- head(delind,1) - 1
### Reorganizing data set for analysis--only vaccinated
dataset <- NULL

for (t in 1:seas2){
  totcases <- sum(studydata[t,2:7])
  totnoncases <- sum(studydata[t,8:13])
  oddscontrol <- totcases/totnoncases * ccratio
  pcontrol <- oddscontrol * totnoncases/Ntot 
  
  for (k in seq_along(vaccdelim[-1])){
    ncases <- studydata[t,2 + k]
    nnoncases <- studydata[t,8 + k]
    
    ncontrols <- rbinom(1,nnoncases,pcontrol)
    # ncontrols <- nnoncases
    if ((ncases > 0 & ncontrols > 0) & (!is.na(ncases) & !is.na(ncontrols))){
      datarec <- c(t,k,1,ncases)
      dataset <- rbind(dataset,datarec,deparse.level = 0)
      datarec <- c(t,k,0,ncontrols)
      dataset <- rbind(dataset,datarec,deparse.level = 0)
    }
  }
}

colnames(dataset) <- c('time','sincevacc','case','count')
dataset <- data.frame(dataset)
dataset$sincevacc <- factor(dataset$sincevacc)

# delind <- which(dataset$case==1 & dataset$count==0)
# dataset <- dataset[-delind,]

cond_logist <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset,method = 'approximate')
summary(cond_logist)

logist <- glm(case ~ sincevacc ,weights = count, data = dataset,family = binomial(link = 'logit'))
summary(logist)
######################################################################################################
######################################################################################################
plot(Ils)
