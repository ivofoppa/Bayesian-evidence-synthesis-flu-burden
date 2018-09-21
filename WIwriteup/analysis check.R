dataset2 <- NULL

for (k in seq(1,length(dataset[,1]),2)){
  cases <- dataset$count[k]
  controls <- dataset$count[k + 1]
  datarec <- c(t,dataset$sincevacc[k],cases,controls)
  dataset2 <- rbind(dataset2,datarec,deparse.level = 0)
}

colnames(dataset2) <- c('time','sincevacc','cases','controls')
dataset2 <- data.frame(dataset2)
dataset2$sincevacc <- factor(dataset2$sincevacc)

logist2 <- glm(cbind(cases,controls) ~ sincevacc, data = dataset2,family = binomial(link = 'logit'))
summary(logist2)

logist <- glm(case ~ sincevacc, data = dataset,weight=count,family = binomial(link = 'logit'))
summary(logist)
