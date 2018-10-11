studydata <- data.frame(studydata0)
studydata1 <- data.frame(studydata0)
studydata12 <- data.frame(studydata10)
studydata2 <- data.frame(studydata20)

delind <- which(rowSums(studydata0[,2:(nvaccat + 1)])==0 | any(is.na(studydata0)))
delind1 <- which(rowSums(studydata0[,2:3]) < 5 | rowSums(studydata0[,4:(nvaccat + 1)]) < 5 | any(is.na(studydata0)))
delind2 <- which(rowSums(studydata20[,2:(nvaccat + 1)])==0 | any(is.na(studydata20)))

if (length(delind)!=0){
  studydata <- studydata[-delind,]
  studydata1 <- studydata1[-delind,]
  studydata12 <- studydata12[-delind,]
  studydata2 <- studydata2[-delind2,]
}
## Defining time-since-vacc delimiters
### Reorganizing data set for analysis--only vaccinated
dataset <- dataset1 <- dataset12 <- dataset2 <- dataset3 <- NULL

seas21 <- studydata$X1

for (t in 1:length(seas21)){
  totcases <- sum(studydata[t,3:(nvaccat + 2)])
  totnoncases <- sum(studydata[t,(nvaccat + 4) : (1 * nvaccat + 3)])
  oddscontrol <- totcases/totnoncases * ccratio
  pcontrol <- ifelse(oddscontrol * totnoncases > Ntot, 1, oddscontrol * totnoncases/Ntot)

  for (k in seq_along(vaccdelim[-1])){
    ncases <- studydata[t,2 + k]
    nnoncases <- studydata[t,nvaccat + 3 + k]
    
    ncontrols <- rbinom(1,nnoncases,pcontrol)

    if ((ncases >= 5 & ncontrols >= 5) & (!is.na(ncases) & !is.na(ncontrols))){
      datarec <- c(seas21[t],k,1,ncases)
      dataset <- rbind(dataset,datarec,deparse.level = 0)
      datarec <- c(seas21[t],k,0,ncontrols)
      dataset <- rbind(dataset,datarec,deparse.level = 0)
    }
  }
  
  for (k in 1:nvaccat2){
    if (k < nvaccat2){
      ncases <- studydata[t,2 + k]
      nnoncases <- studydata[t,nvaccat + 3 + k]

      ncontrols <- rbinom(1,nnoncases,pcontrol)
    } else {
      ncases <- sum(studydata[t,((nvaccat2 + 2) : (nvaccat + 2))])
      nnoncases <- sum(studydata[t,(nvaccat + nvaccat2 + 3):((2 * nvaccat + 3))])
      
      ncontrols <- rbinom(1,nnoncases,pcontrol)
    } 
    
    if ((ncases >= ncrit & ncontrols >= ncrit) & (!is.na(ncases) & !is.na(ncontrols))){
      datarec <- c(seas21[t],k,1,ncases)
      dataset1 <- rbind(dataset1,datarec,deparse.level = 0)
      datarec <- c(seas21[t],k,0,ncontrols)
      dataset1 <- rbind(dataset1,datarec,deparse.level = 0)
    }
  }
}

seas22 <- studydata2$X1

for (t in 1:length(seas22)){
  totcases2 <- sum(studydata2[t,3:(nvaccat + 2)])
  totnoncases2 <- sum(studydata2[t,(nvaccat + 4) : (1 * nvaccat + 3)])
  oddscontrol2 <- totcases2/totnoncases2 * ccratio
  pcontrol2 <- ifelse(oddscontrol2 * totnoncases2 > Ntot, 1, oddscontrol2 * totnoncases2/Ntot)
  
  for (k in seq_along(vaccdelim[-1])){
    ncases2 <- studydata2[t,2 + k]
    nnoncases2 <- studydata2[t,nvaccat + 3 + k]
    
    ncontrols2 <- rbinom(1,nnoncases2,pcontrol2)
    
    if ((ncases2 >= 5 & ncontrols2 >= 5) & (!is.na(ncases2) & !is.na(ncontrols2))){
      datarec2 <- c(seas22[t],k,1,ncases2)
      dataset2 <- rbind(dataset2,datarec2,deparse.level = 0)
      datarec2 <- c(seas22[t],k,0,ncontrols2)
      dataset2 <- rbind(dataset2,datarec2,deparse.level = 0)
    }
  }
  
  for (k in 1:nvaccat2){
    if (k < nvaccat2){
      ncases2 <- studydata2[t,2 + k]
      nnoncases2 <- studydata2[t,nvaccat + 3 + k]
      
      ncontrols2 <- rbinom(1,nnoncases2,pcontrol2)
    } else {
      ncases2 <- sum(studydata2[t,((nvaccat2 + 2) : (nvaccat + 2))])
      nnoncases2 <- sum(studydata2[t,(nvaccat + nvaccat2 + 3):((2 * nvaccat + 3))])
      
      ncontrols2 <- rbinom(1,nnoncases2,pcontrol2)
    } 

    if ((ncases2 >= ncrit & ncontrols2 >= ncrit) & (!is.na(ncases2) & !is.na(ncontrols2))){
      datarec2 <- c(seas22[t],k,1,ncases2)
      dataset12 <- rbind(dataset12,datarec2,deparse.level = 0)
      datarec2 <- c(seas22[t],k,0,ncontrols2)
      dataset12 <- rbind(dataset12,datarec2,deparse.level = 0)
    }
  }
}

dataset3 <- dataset
colnames(dataset3) <- c('time','sincevacc','case','count')
dataset3 <- data.frame(dataset3)

colnames(dataset) <- c('time','sincevacc','case','count')
dataset <- data.frame(dataset)
sincevacc <- dataset$sincevacc
dataset$sincevacc <- factor(dataset$sincevacc)

colnames(dataset1) <- c('time','sincevacc','case','count')
dataset1 <- data.frame(dataset1)
sincevacc <- dataset1$sincevacc
dataset1$sincevacc <- factor(dataset1$sincevacc)

colnames(dataset12) <- c('time','sincevacc','case','count')
dataset12 <- data.frame(dataset12)
sincevacc <- dataset12$sincevacc
dataset12$sincevacc <- factor(dataset12$sincevacc)

colnames(dataset2) <- c('time','sincevacc','case','count')
dataset2 <- data.frame(dataset2)
dataset2$sincevacc <- factor(dataset2$sincevacc)
######################################################################################################
##  Crude analysis, only adjusting for time ##########################################################
##  Create new data set first               ##########################################################
######################################################################################################
dataset4 <- NULL
VEls <- trueVEls <- NULL
for (t in 1:seas2){
  totcases <- sum(studydata[t,2:(nvaccat + 2)])
  totnoncases <- sum(studydata[t,(nvaccat + 3) : (2 * nvaccat + 3)])
  oddscontrol <- totcases/totnoncases * ccratio
  pcontrol <- min(oddscontrol * totnoncases/Ntot,1)
  
  ncasesnv <- studydata[t,2]
  nnoncasesnv <- studydata[t,(nvaccat + 3)]
  nnoncasesnv2 <- studydata12[t,(nvaccat + 3)]
  
  ncasesv <- sum(studydata[t,3:(nvaccat + 2)])
  nnoncasesv <- sum(studydata[t,(nvaccat + 4) : (1 * nvaccat + 3)])
  nnoncasesv2 <- sum(studydata12[t,(nvaccat + 4) : (1 * nvaccat + 3)])
  
  ncontrolsv <- rbinom(1,nnoncasesv,pcontrol)
  ncontrolsnv <- rbinom(1,nnoncasesnv,pcontrol)
  ncontrolsv2 <- nnoncasesv2
  ncontrolsnv2 <- nnoncasesnv2

  late <- ifelse(t > seas2/2,1,0)
  if (ncasesv >= 5 & ncasesnv >= 5){
    datarec <- c(t,1,ncasesv,ncontrolsv,late)
    dataset4 <- rbind(dataset4,datarec,deparse.level = 0)
    datarec <- c(t,0,ncasesnv,ncontrolsnv,late)
    dataset4 <- rbind(dataset4,datarec,deparse.level = 0)
    
    veest <- 1 - ncasesv * ncontrolsnv / (ncasesnv * ncontrolsv)
    VEls <- c(VEls, veest)
    trueveest <- 1 - ncasesv * ncontrolsnv2 / (ncasesnv * ncontrolsv2)
    trueVEls <- c(trueVEls, trueveest)
  }
}
colnames(dataset4) <- c('time','vacc','cases','controls','late')
dataset4 <- data.frame(dataset4)
