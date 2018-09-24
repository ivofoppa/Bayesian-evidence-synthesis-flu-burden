

ls1 <- sapply(1:5, function(x) studydata[tm,2 + x])
ls2 <- sapply(1:5, function(x) studydata[tm,8 + x])

oddsls <- ls1/ls2
orls <- oddsls/oddsls[1]

tm <- 20
dta <- dataset[which(dataset$time == tm),]

plot(rowSums(studydata[,2:7]),type = 'l')
