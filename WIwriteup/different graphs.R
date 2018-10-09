filename <- paste0('data_National_2010_16_5ag.dat')
setwd('..')
setwd("./RSVdata")

datarr <- read.table(file = filename, header = T)
datarr5 <- datarr[which(datarr$age==5),]

rcu <- datarr5$rcu
RSV <- datarr5$RSV
rsvstd <- RSV/(max(RSV))*max(rcu - min(rcu)) + min(rcu)

H3 <- datarr5$AH3
H3std <- H3/(max(H3))*max(rcu - min(rcu)) + min(rcu)

corrdata <- data.frame(rcu,rsvstd,H3std)

corrdata <- corrdata[order(corrdata$rcu),]

plot(corrdata$rcu,corrdata$rsvstd)
plot(corrdata$rcu,corrdata$H3std)

###################################################################################################
variables7 <- c(paste0('muB0[3:',N,']'),'b20','b21','b22')
###################################################################################################

codaarr <- data.frame(codaarr)

mu <- sapply(3:N, function(x) eval(parse(test = paste0('median(codaarr$muB0.',x,'.)'))))

b20 <- median(codaarr$b20)
b21 <- median(codaarr$b21)
b22 <- median(codaarr$b22)

rsv2 <- RSV[1:(N-2)]
rsv1 <- RSV[2:(N-1)]
rsv0 <- RSV[3:N]

rsvEM <- (rsv0*b20 + rsv1*b21 + rsv2*b22) * pop[3:N]

plot(mort[-c(1:2)])
lines(mu + rsvEM,col = 'red')

## constructing alternative baselines
mu1 <- rcu[3:N] - .1*rsvEM
mu2 <- rcu[3:N] - .5*rsvEM
mu3 <- rcu[3:N] - rsvEM
mu4 <- rcu[3:N] - 2*rsvEM
mu5 <- rcu[3:N] - 5*rsvEM
mu6 <- rcu[3:N] - 10*rsvEM

cols <- rainbow(7)
setwd('../') ## change to root directory
setwd('./RSVwriteup')
pdf('BL.pdf',paper='USr',height = 8.5,width = 11) 
plot(mu1,type = 'l',ylim = c(min(c(mu1,mu2,mu3,mu4,mu5,mu6)),max(c(mu1,mu2,mu3,mu4,mu5,mu6))),col = cols[1], ylab = '',xlab = 'Index Week')

for (k in 2:6){
  eval(parse(text = paste0('lines(mu',k,',col = cols[',k,'])')))
}

rsv <- RSV/max(RSV)*max(c(mu1,mu2,mu3,mu4,mu5,mu6))*.3

lines(rsv[3:N],col = cols[7])

legend('bottom',c('BL fact 0.1','BL fact 0.5','BL fact 1','BL fact 2','BL fact 5','BL fact 10','RSV'), lty = 1, col = cols,bty = 'n',
       ncol = 5)
dev.off()
