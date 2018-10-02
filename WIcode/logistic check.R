dataset3 <- NULL

early_crit <- 75

for (t in 1:seas2){
  
  lateseas <- ifelse(t < early_crit, 0,1)  
  for (vacc in 0:1){
    if (vacc==0){
      cases <- studydata[t,2]
      controls <- studydata[t,8]
    } else {
      cases <- sum(studydata[t,3:7])
      controls <- sum(studydata[t,8:13])
    }
    
    if ((cases > 0 & controls > 0) & (!is.na(cases) & !is.na(controls))){
      datarec <- c(t,vacc, lateseas, cases,controls)
      dataset3 <- rbind(dataset3,datarec,deparse.level = 0)
    }
  }
}
colnames(dataset3) <- c('time','vacc','lateseas','cases','controls')
dataset3 <- data.frame(dataset3)

cases <- dataset3$cases
controls <- dataset3$controls
vacc <- dataset3$vacc
lateseas <- dataset3$lateseas

logitmod <- glm(cbind(cases,controls) ~ vacc + lateseas + lateseas*vacc, family = binomial(link = 'logit'))

summary(logitmod)
