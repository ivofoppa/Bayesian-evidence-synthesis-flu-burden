library(stats)
library(RColorBrewer)

bfolder <- "C:/Users/vor1/Dropbox/Misc work/Waning immunity/WI git project/"
setwd(paste0(bfolder,"WIgraphs"))
### Leaky: run leaky code first ('waning 1 strain leaky VE range.R')
#########################################################################################
## Plot the VE estimates by VE  #########################################################
#########################################################################################
### Individual VE=.5
k <- 5
VEls <- VEarr[k,]
casels <- casearr[k,]

### Without cases
png(filename = 'Leaky_VEest_50.png')
selind <- which(casearr[1,] > 0)
plot(head(times,-1),sapply(1:length(VEls),function(x) ifelse(x %in% selind, VEls[x],NA)),type = 'l',
     ylim = c(min(0,min(VEls)),max(VEls)*1.1),ylab = 'VE',xlab = 'Day',
     main = 'VE Estimates by True VE',lwd = 1,col = 'darkorchid')

lines(head(times,-1),rep(k/10,length(times)-1),col = 'blue')
lines(head(times,-1),sapply(1:length(VEls),function(x) ifelse(x %in% selind, VEls[x],NA)),col = 'darkorchid')

legend('left',c('VE est.','True VE','         '),col = c('darkorchid','blue','white'), 
       lty = 1,bty = 'n',lwd = 1,seg.len = 2.5)
dev.off()
### With cases
png(filename = 'Leaky_VEest_50_cases.png')
selind <- which(casearr[1,] > 0)
plot(head(times,-1),sapply(1:length(VEls),function(x) ifelse(x %in% selind, VEls[x],NA)),type = 'l',
     ylim = c(min(0,min(VEls)),max(VEls)*1.1),ylab = 'VE',xlab = 'Day',
     main = 'VE Estimates by True VE',lwd = 1,col = 'darkorchid')

lines(head(times,-1),rep(k/10,length(times)-1),col = 'blue')
lines(head(times,-1),sapply(1:length(VEls),function(x) ifelse(x %in% selind, VEls[x],NA)),col = 'darkorchid')
lines(head(times,-1),casels/max(casels)*max(VEls)/2,col = 'red')

legend('left',c('VE est.','True VE','Incidence'),col = c('darkorchid','blue','red'), 
       lty = 1,bty = 'n',lwd = 1,seg.len = 2.5)
dev.off()
#########################################################################################
#########################################################################################
### Range of VE
k <- 1
VEls <- VEarr[k,]
#cols <- brewer.pal(3, 'Spectral')
cols <- c('red','green','blue')
png(filename = 'Leaky_VEest_range.png')
selind <- which(casearr[1,] > 0)
plot(head(times,-1),sapply(1:length(VEls),function(x) ifelse(x %in% selind, VEls[x],NA)),type = 'l',
     ylim = c(min(0,min(VEls)),max(VEarr)*1.3),ylab = 'VE est.',xlab = 'Day',
     main = 'VE Estimates by True VE',lwd = 1,col = cols[1])

for(k in 2:3){
  selind <- which(casearr[k,] > 0)
  lines(head(times,-1),sapply(1:length(VEls),function(x) ifelse(x %in% selind, VEarr[k,x],NA)),col = cols[k])
}

VElabs <- sapply(VErnge, function(x) paste0(x*100))

legend('top',VElabs,col = cols, title = 'True VE (%)',horiz = T,
       lty = 1,bty = 'n',lwd = 1,seg.len = .5)
dev.off()
#########################################################################################
#########################################################################################
## Plot the epi curves by VE  ###########################################################
#########################################################################################
k <- 1
casels <- casearr[k,]

cols <- brewer.pal(9, 'Spectral')

png(filename = 'Epicurve.png')
plot(head(times,-1),casels*100,type = 'l',ylim = c(min(0,min(casels)),max(casels*100)*1.1),ylab = 'Inc. dens.',xlab = 'Day',
     main = 'Epi Curves by True VE',lwd = 1,col = cols[1])

for(k in 2:9){
  lines(head(times,-1),casearr[k,]*100,col = cols[k])
}

VElabs <- sapply(1:9, function(x) paste0(x*10))

legend('left',VElabs,col = cols, title = 'VE (%)',
       lty = 1,bty = 'n',lwd = 1,seg.len = 2.5)
dev.off()
#########################################################################################
##  Two viruses with all successfully vacc. immune to virus 2 ###########################
## ('waning 2 strain full alpha 1.R')                         ###########################
##  R01 = 1.25, R02 = 1.75, vaccination uptake = 0.3          ###########################
#########################################################################################
casesv <- -diff(xv1 + xv2)
casesnv <- -diff(xnv)

selind <- which((casesv + casesnv)>0)

infprsv <- (yv11 + yv21 + ynv1)*beta1 + (1 - xv1/(xv1 + xv2))*(yv22 + ynv2)*beta2
infprsnv <- (yv11 + yv21 + ynv1)*beta1 + (yv22 + ynv2)*beta2
rls <- casesv/casesnv
## odds of vacc among cases
VEls <- 1-( rls/(vacc/(1-vacc)))

trueVEls <- 1-(infprsv/infprsnv)[-1]
Nsusc <- (xv1 + xv2 + xnv)[-1]
VEmax <- max(VEls,trueVEls)

a3 <- ceiling(max(Nsusc)) 

y3tck2 <- seq(0,a3,a3/5)/a3*VEmax
y3tck2lab <- sapply(round(seq(0,a3,a3/5)/1000000)*1000000, function(x) toString(x))

times1 <- times[-1]
par(mar = c(5,4,4,5))
png(filename = 'two_viruses1.png')
plot(times1,VEls,type = 'l',ylim = c(min(VEls,0),VEmax),ylab = 'VE',xlab = 'Day',
     main = '',lwd = 2)
lines(times1,trueVEls, col = 'white')
lines(times1,Nsusc/a3*VEmax, col = 'white',lwd = 2)
axis(side = 4,at= y3tck2,labels = y3tck2lab)
mtext('Susceptibles', 4, line = 3)

legend('bottomleft',c('VE','       ','          '),col = c('black','white','white'),
       lty = 1,bty = 'n',lwd = 2,seg.len = 1)
dev.off()
#########################################################################################
times1 <- times[-1]
png(filename = 'two_viruses2.png')
plot(times1,VEls,type = 'l',ylim = c(min(VEls,0),VEmax),ylab = 'VE',xlab = 'Day',
     main = '',lwd = 2)
lines(times1,trueVEls, col = 'red')
lines(times1,Nsusc/a3*VEmax, col = 'white',lwd = 2)
axis(side = 4,at= y3tck2,labels = y3tck2lab)
mtext('Susceptibles', 4, line = 3)

legend('bottomleft',c('VE','True VE','          '),col = c('black','red','white'),
       lty = 1,bty = 'n',lwd = 2,seg.len = 1)
dev.off()
#########################################################################################
png(filename = 'two_viruses3.png')
plot(times1,VEls,type = 'l',ylim = c(min(VEls,0),VEmax),ylab = 'VE',xlab = 'Day',
     main = '',lwd = 2)
lines(times1,trueVEls, col = 'red')
lines(times1,Nsusc/a3*VEmax, col = 'blue',lwd = 2)
axis(side = 4,at= y3tck2,labels = y3tck2lab)
mtext('Susceptibles', 4, line = 3)

legend('bottomleft',c('VE','True VE','Pop. susc.'),col = c('black','red','blue'),
       lty = 1,bty = 'n',lwd = 2,seg.len = 1)
dev.off()
#########################################################################################
###  Virus distribution #################################################################
#########################################################################################
png(filename = 'prop_virus1.png')

plot(times,(yv11 + yv21 + ynv1)/(yv11 + yv21 + ynv1 + yv22 + ynv2),type = 'l' ,ylim = c(0.4,1),
     main = '',
     ylab = 'Virus 1',xlab = 'Day')
dev.off()

png(filename = 'prop_virus1b.png')
propvir1 <- (yv11 + yv21 + ynv1)/(yv11 + yv21 + ynv1 + yv22 + ynv2)
propsuscv2 <- xv2/(xv1 + xv2)
plot(times,propvir1,type = 'l' ,ylim = c(0,1), lwd = 2,
     main = 'Proportion of Virus 1',
     ylab = 'Virus 1',xlab = 'Day')
lines(propsuscv2/max(propsuscv2)*max(propvir1),col = 'green',lwd = 2)
dev.off()
