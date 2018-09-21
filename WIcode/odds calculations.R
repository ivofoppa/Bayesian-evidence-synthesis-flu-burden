ls1 <- sapply(1:5, function(x) sum(studydata[,1 + x]))
ls2 <- sapply(1:5, function(x) sum(studydata[,7 + x]))

oddsls <- ls1/ls2
orls <- oddsls/oddsls[1]
 