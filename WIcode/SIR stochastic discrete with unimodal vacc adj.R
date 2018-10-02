## Discrete stochastic transmission model with unimodal vaccination uptake and adjustment
library(binhf)
library(survival)
Ntot <- 2000000
r0 <- 1.6
delta <- 1/4 ## infectious period
beta <- r0 * delta

phi <- .6

pSE <- function(r0,I,phi,vacc){ ## infection probability
  lambda <- beta * (1-phi)^vacc * I / Ntot
  return (as.numeric(1-exp(-lambda)))
}

seas <- 500 ## Number of days in epidemic

vacc1 <- 0.47

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
Invinit <- 0
Snvinit <- Ntot
Rnvinit <- 0

Ivinit <- rep(0,seas + prevdur)
Svinit <- rep(0,seas + prevdur)
Rvinit <- rep(0,seas + prevdur)

inits <- list(Svinit=Svinit,Snvinit=Snvinit,
              Ivinit=Ivinit,Invinit=Invinit,Rvinit=Rvinit,Rnvinit=Rnvinit)

# sim <- 1
ccratio <- 3

Snv <- as.integer(inits$Snvinit)
Inv <- as.integer(inits$Invinit)
Rnv <- as.integer(inits$Rnvinit)

Sv <- as.integer(inits$Svinit)
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
studydata2 <- array(0,dim = c(seas,1 + 2 * nvaccat + 2))

infnum <- 100

pinfnv <-  Snv/(Snv + (1 - phi)*sum(Sv))  ### probability that infections are non-vaccinated
infnv <- rbinom(1,infnum,pinfnv)
Inv <- Inv + infnv

infv <- infnum - infnv
pSv <- Sv/max(sum(Sv),1)
Iv <- Iv + c(rmultinom(1,infv,pSv))

I <- sum(c(Iv,Inv))

Ils <- NULL

while (time <= seas + prevdur){
  pnvinf <- pSE(r0,I,phi,0)
  pvinf <- pSE(r0,I,phi,1)
  
  newinfnv <- rbinom(1, Snv, pnvinf)
  newinfv <- sapply(Sv, function(sv) rbinom(1, sv, pvinf))
  
  newinfv2 <- sapply(seq_along(newinfv), function(k) ifelse(k <= (time - prevdur), 0, newinfv[k]))
  
  newremnv <- rbinom(1,Inv,1 - exp(-delta))
  newremv <- sapply(Iv, function(iv) rbinom(1,iv, 1 - exp(-delta)))

  Ils <- c(Ils,sum(newinfv) + newinfnv)  
  Snv <- Snv - newinfnv
  Sv <- Sv - newinfv
  
  Inv <- Inv + newinfnv - newremnv
  Iv <- Iv + newinfv - newremv
  
  I <- sum(c(Iv,Inv))
  
  Rnv <- Rnv + newremnv
  Rv <- Rv + newremv
  
  ###################################################################################################
  ### "Study" is conducted ##########################################################################
  ###################################################################################################
  casesls <- c(newinfnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(newinfv[vaccdelim[x] : vaccdelim[x + 1]])))
  casesls2 <- c(newinfnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(newinfv2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  controlsvsum <- Sv + Rv + Iv
  controlsvsum2 <- sapply(seq_along(controlsvsum), function(k) ifelse(k <= (time - prevdur), 0, controlsvsum[k]))

  controlsls <- c(Snv + Rnv + Inv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum[vaccdelim[x] : vaccdelim[x + 1]])))
  controlsls2 <- c(Snv + Rnv + Inv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  studydata[time - prevdur,] <- c(time - prevdur,casesls,controlsls)
  studydata2[time - prevdur,] <- c(time - prevdur,casesls2,controlsls2)
  ###################################################################################################
  ### Stop when no more transmission 
  if (sum(Iv)==0 & Inv == 0){
    break()
  }
  ###################################################################################################
  ### Vaccination: proceeding time since vacc. and adding new vaccinees
  Sv <- shift(Sv,1)
  svacc <- rbinom(1,Snv,1 - exp(-v[time]))
  Sv[1] <- svacc; Snv <- Snv - svacc
  
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
studydata2 <- studydata2[-delind,]

seas2 <- head(delind,1) - 1
### Reorganizing data set for analysis--only vaccinated
dataset <- dataset2 <- NULL

for (t in 1:seas2){
  totcases <- sum(studydata[t,2:7])
  totnoncases <- sum(studydata[t,8:13])
  oddscontrol <- totcases/totnoncases * ccratio
  pcontrol <- oddscontrol * totnoncases/Ntot 
  
  totcases2 <- sum(studydata2[t,2:7])
  totnoncases2 <- sum(studydata2[t,8:13])
  oddscontrol2 <- totcases2/totnoncases2 * ccratio
  pcontrol2 <- oddscontrol2 * totnoncases2/Ntot 
  
  for (k in seq_along(vaccdelim[-1])){
    ncases <- studydata[t,2 + k]
    nnoncases <- studydata[t,8 + k]
    
    ncases2 <- studydata2[t,2 + k]
    nnoncases2 <- studydata2[t,8 + k]
    
    ncontrols <- rbinom(1,nnoncases,pcontrol)
    ncontrols2 <- rbinom(1,nnoncases2,pcontrol2)
    # ncontrols <- nnoncases
    if ((ncases > 0 & ncontrols > 0) & (!is.na(ncases) & !is.na(ncontrols))){
      datarec <- c(t,k,1,ncases)
      dataset <- rbind(dataset,datarec,deparse.level = 0)
      datarec <- c(t,k,0,ncontrols)
      dataset <- rbind(dataset,datarec,deparse.level = 0)
    }
    if ((ncases2 > 0 & ncontrols2 > 0) & (!is.na(ncases2) & !is.na(ncontrols2))){
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

filepath <- paste0('C:/Users/VOR1/Documents/GitHub/Waning-Immunity-artefact/WIwriteup/WIplots/simul_adj_',prevdur,'_',vdur,'.RData')
save(dataset,dataset2,studydata,studydata2,file = filepath)

logist <- glm(case ~ sincevacc ,weights = count, data = dataset,family = binomial(link = 'logit'))
summary(logist)
######################################################################################################
######################################################################################################
plot(Ils)
