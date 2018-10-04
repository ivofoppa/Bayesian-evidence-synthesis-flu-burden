studydata <- studydata0
studydata2 <- studydata20

delind <- which(rowSums(studydata0[,2:(nvaccat + 1)])==0 | any(is.na(studydata0)))
delind2 <- which(rowSums(studydata20[,2:(nvaccat + 1)])==0 | any(is.na(studydata20)))

studydata <- studydata0[-delind,]
studydata2 <- studydata20[-delind2,]

seas2 <- head(delind,1) - 1
## Defining time-since-vacc delimiters
### Reorganizing data set for analysis--only vaccinated
dataset <- dataset2 <- NULL

for (t in 1:seas2){
  totcases <- sum(studydata[t,2:(nvaccat + 2)])
  totnoncases <- sum(studydata[t,(nvaccat + 3) : (1 * nvaccat + 3)])
  oddscontrol <- totcases/totnoncases * ccratio
  pcontrol <- ifelse(oddscontrol * totnoncases > Ntot, 1, oddscontrol * totnoncases/Ntot)
  
  totcases2 <- sum(studydata2[t,2:(nvaccat + 2)])
  totnoncases2 <- sum(studydata2[t,(nvaccat + 3) : (1 * nvaccat + 3)])
  oddscontrol2 <- totcases2/totnoncases2 * ccratio
  pcontrol2 <- ifelse(oddscontrol2 * totnoncases2 > Ntot, 1, oddscontrol2 * totnoncases2/Ntot)
  
  for (k in seq_along(vaccdelim[-1])){
    ncases <- studydata[t,2 + k]
    nnoncases <- studydata[t,nvaccat + 2 + k]
    
    ncases2 <- studydata2[t,2 + k]
    nnoncases2 <- studydata2[t,nvaccat + 2 + k]
    
    ncontrols <- rbinom(1,nnoncases,pcontrol)
    ncontrols2 <- rbinom(1,nnoncases2,pcontrol2)
    
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
######################################################################################################
##  Crude analysis, only adjusting for time ##########################################################
##  Create new data set first               ##########################################################
######################################################################################################
dataset3 <- NULL
VEls <- NULL
ccratio <- 2
for (t in 1:seas2){
  totcases <- sum(studydata[t,2:(nvaccat + 2)])
  totnoncases <- sum(studydata[t,(nvaccat + 3) : (1 * nvaccat + 3)])
  oddscontrol <- totcases/totnoncases * ccratio
  pcontrol <- oddscontrol * totnoncases/Ntot
  
  ncasesnv <- studydata[t,2]
  nnoncasesnv <- studydata[t,(nvaccat + 3)]
  
  ncasesv <- sum(studydata[t,3:(nvaccat + 2)])
  nnoncasesv <- sum(studydata[t,(nvaccat + 4) : (1 * nvaccat + 3)])
  
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
