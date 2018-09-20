## Discrete stochastic transmission model with contact matrix
## Latent period (from Comflu_bas): 30% 1 day; 50% 2 days and 20% 3 days;
## Infectious period: 30% 3, 40% 4, 20% 5, 10% 6 days
library(binhf)
library(survival)
Ntot <- 2000000
## Creating population according to latent (rows) and infectious periods (columns) in each of the age groups
r0 <- 1.5
gamma <- 1/2 ## latent period
delta <- 1/4 ## infectious period
beta <- r0 * delta

phi <- .6

pSE <- function(r0,I,phi,vacc){ ## infection probability
  lambda <- beta * (1-phi)^vacc * I / Ntot
  return (as.numeric(1-exp(-lambda)))
}

seas <- 150 ## Number of days in epidemic
### Generating the vaccination rate over time
vacc <- .4

campdur <- 50
vrate <- -log(1 - vacc)/campdur

v <- c(.10,rep(vrate,campdur),rep(0,seas - 51))
vaccdelim <- c(seq(1,60,12),seas)
nvaccat <- length(vaccdelim) ## number of time-since-vacc categories

Envinit <- 10
Invinit <- 0
Snvinit <- Ntot - Envinit - Invinit
Rnvinit <- 0

Evinit <- rep(0,seas)
Ivinit <- rep(0,seas)
Svinit <- rep(0,seas)
Rvinit <- rep(0,seas)

inits <- list(Svinit=Svinit,Snvinit=Snvinit,Evinit=Evinit,Envinit=Envinit,
              Ivinit=Ivinit,Invinit=Invinit,Rvinit=Rvinit,Rnvinit=Rnvinit)

nsim <- 10

sim <- 1

while ( sim <= nsim ){

  Snv <- as.integer(inits$Snvinit)
  Env <- as.integer(inits$Envinit)
  Inv <- as.integer(inits$Invinit)
  Rnv <- as.integer(inits$Rnvinit)
  
  Sv <- as.integer(inits$Svinit)
  Ev <- as.integer(inits$Evinit)
  Iv <- as.integer(inits$Ivinit)
  Rv <- as.integer(inits$Rvinit)
  
  I <- sum(Iv + Inv)
  
  Snvls <- Svls <- Envls <- Evls <- Invls <- Ivls <- Rnvls <- Rvls <- NULL
  
  time <- 1
  
  studydata <- array(0,dim = c(seas,nvaccat + nvaccat + 1))

  while (time <= seas){
    pnvinf <- pSE(r0,I,phi,0)
    pvinf <- pSE(r0,I,phi,1)
    
    newlatnv <- rbinom(1, Snv, pnvinf)
    newlatv <- sapply(Sv, function(sv) rbinom(1, sv, pvinf))

    newinfnv <- rbinom(1,Env,1 - exp(-gamma))
    newinfv <- sapply(Ev,function(ev) rbinom(1,ev,1 - exp(-gamma)))

    newremnv <- rbinom(1,Inv,1 - exp(-delta))
    newremv <- sapply(Iv, function(iv) rbinom(1,iv, 1 - exp(-delta)))

    Snv <- Snv - newlatnv
    Sv <- Sv - newlatv
    
    Env <- Env + newlatnv - newinfnv
    Ev <- Ev + newlatv - newinfv

    Inv <- Inv + newinfnv - newremnv
    Iv <- Iv + newinfv - newremv

    I <- sum(Iv + Inv)
    
    Rnv <- Rnv + newremnv
    Rv <- Rv + newremv
    
    Invls <- c(Invls,Inv)
    Ivls <- c(Ivls,sum(Iv))
    ###################################################################################################
    ### "Study" is conducted ##########################################################################
    ###################################################################################################
    casesvls <- sapply(seq_along(vaccdelim[-1]), function(x) sum(newinfv[vaccdelim[x] : vaccdelim[x + 1]]))
    
    controlsvsum <- Sv + Ev + Rv
    controlsvls <- sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum[vaccdelim[x] : vaccdelim[x + 1]]))
    
    studydata[time,] <- c(time,casesvls,newinfnv,controlsvls, Snv + Env + Rnv)
    ###################################################################################################
    ###################################################################################################
    ### Vaccination
    Sv <- shift(Sv,1)
    svacc <- rbinom(1,Snv,1 - exp(-v[time]))
    Sv[1] <- svacc; Snv <- Snv - svacc
    
    Ev <- shift(Ev,1)
    evacc <- rbinom(1,Env,1 - exp(-v[time]))
    Ev[1] <- evacc; Env <- Env - evacc
    
    Iv <- shift(Iv,1)
    ivacc <- rbinom(1,Inv,1 - exp(-v[time]))
    Iv[1] <- ivacc; Inv <- Inv - ivacc

    time <- time + 1
  }
  
  ### Reorganizing data set for analysis--only vaccinated
  dataset <- NULL
  
  for (t in 1:seas){
    for (k in seq_along(vaccdelim[-1])){
      datarec <- c(t,k,1,studydata[t,1 + k])
      dataset <- rbind(dataset,datarec,deparse.level = 0)
      datarec <- c(t,k,0,studydata[t,7 + k])
      dataset <- rbind(dataset,datarec,deparse.level = 0)
    }
  }
  
  colnames(dataset) <- c('time','sincevacc','case','count')
  dataset <- data.frame(dataset)
  dataset$sincevacc <- factor(dataset$sincevacc)
  
  delind <- which(dataset$count==0)
  dataset <- dataset[-delind,]
  
  cond_logist <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset,method = 'approximate')
  
  assign(paste0('Infectious',sim),Ils)
  sim <- sim + 1
}

cols <- rainbow(nsim)
xlim <- 110
xlow <- 30
plot(Infectious1[xlow:xlim], type = 'l',col = cols[1],xlab = 'Outbreak Day',ylab = 'Number Infectious')

#plot(Ils[1:150], type = 'l',col = cols[1])

for (k in 2:nsim){
  eval(parse(text = paste0('lines(infectious',k,'[xlow:xlim],col = cols[k], type = \"l\")')))
}
######################################################################################################
######################################################################################################
