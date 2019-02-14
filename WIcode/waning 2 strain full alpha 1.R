### Investigating the presence of waning when all-or none effect with against specific influenza strains
library(deSolve)
bfolder <- 'C:/Users/VOR1/Documents/GitHub/' ## Define root path (where repository)
## asign values to parameters
Ntot <- 20000000
vacc <- 0.3
# phi <- 1
Nv0 <- round(Ntot*vacc)
Nnv0 <- Ntot - Nv0
alpha <- 1 #VE

gamma <- 0.333

beta1 <- 1.25*gamma/(Ntot)
beta2 <- 1.75*gamma/(Ntot)
### Initial conditions
yv10 <- yv1 <- 1; 

ynv10 <- ynv1 <- 1; 
ynv20 <- ynv2 <- 1; 

xv0 <- Nv0 - yv10
xnv0 <- Nnv0 - ynv10 - ynv20
z0 <- 0
parameters <- c(gamma=gamma,beta1=beta1, beta2=beta2)
### Define model
KKmod <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dxv <- - beta1*xv*(yv1 + ynv1) 
    dxnv <- - beta1*xnv*(yv1 + ynv1) - beta2*xnv*ynv2
    
    dyv1 <- beta1*xv*(yv1 + ynv1) - yv1*gamma
    
    dynv1 <- beta1*xnv*(yv1 + ynv1) - ynv1*gamma
    dynv2 <- beta2*xnv*ynv2 - ynv2*gamma
    
    dyv1new <- beta1*xv*(yv1 + ynv1)
    
    dynv1new <- beta1*xnv*(yv1 + ynv1)
    dynv2new <- beta2*xnv*ynv2
    
    dz <- (ynv1 + ynv2 + yv1)*gamma
    
    # return the rate of change
    list(c(dxv,dxnv,dynv1,dynv2,dyv1,
           dynv1new,dynv2new,dyv1new,dz))
  })
}
### define vector with initial conditions for ode function
state <- c(xv=xv0,xnv=xnv0,ynv1=ynv10,ynv2=ynv20,yv1=yv10,
           ynv1new=0,ynv2new=0,yv1new=0,z=z0)

seas <- 300
times <- seq(0, seas, by = 0.01)

out <- data.frame(ode(y = state, times = times, func = KKmod, parms = parameters))

xv <- out$xv
xnv <- out$xnv
yv1new <- out$yv1new
ynv1new <- out$ynv1new
ynv2new <- out$ynv2new
###################################################################################################
###################################################################################################
dxv <- -diff(xv)
dxnv <- -diff(xnv)

dynv1 <- diff(ynv1new)
dynv2 <- diff(ynv2new)
dyv1 <- diff(yv1new)

selind <- which((dxv + dxnv)>0)

cases1 <- dyv1 + dynv1; cases1 <- cases1[selind]
cases2 <- dynv2; cases2 <- cases2[selind]

newcases <- cases1 + cases2;
Nsusc <- Ntot - (yv1new + ynv1new + ynv2new)[selind]
casesv <- dyv1; casesv <- casesv[selind]
casesnv <- dynv1 + dynv2; casesnv <- casesnv[selind]

infprsv <- (dyv1 + dynv1)*beta1 + (1-alpha)*(dynv2)*beta2
infprsnv <- (dyv1 + dynv1)*beta1 + (dynv2)*beta2
rls <- casesv/casesnv
## odds of vacc among cases
VEls <- 1-( rls/(vacc/(1-vacc)))

trueVEls <- 1-infprsv/infprsnv

a2 <- ceiling(max(cases2)) 
b2 <- nchar(a2) ## how many digits in number?

a1 <- ceiling(max(cases1)) 
b1 <- nchar(a1) ## how many digits in number?

y2uplim0 <- ceiling(a2/10^(b2-1))
y2uplim <- y2uplim0/10

y2tck2 <- seq(0,y2uplim,y2uplim/5)*10^b2/a2*a1
y2tck2lab <- sapply(round(seq(0,y2uplim,y2uplim/5)*10^b2), function(x) toString(x))
y2toax <- y2uplim0*10^(b2-1)

y1tck2 <- seq(0,y1uplim,.2)*10^b1/a1*a2
y1tck2lab <- sapply(round(seq(0,y1uplim,.2)*10^b1), function(x) toString(x))
y1toax <- y1uplim0*10^(b1-1)

y1uplim0 <- ceiling(a1/(10^(b1-1)))
y1uplim <- ifelse(y1uplim0%%2==0,y1uplim0,y1uplim0+1)/10
y1toax <- y1uplim0*10^(b1-1)

par(mar = c(5,4,4,5))
times1 <- times[selind]
### plot comparing viruses
plot(times1,cases1,type = 'l',ylab = '# cases',xlab = 'Days',col = 'blue',lwd = 2,ylim = c(0,max(cases1)))
# lines(times1,cases1/a1*a2,col = 'red',lwd = 2)
lines(times1,cases2,col = 'red',lwd = 2)

# axis(side = 4,at= y1tck2,labels = y1tck2lab)
# mtext('# cases, virus 2', 4, line = 3)
legend('topleft',c('Virus 1','Virus 2'),col = c('blue','red'),
       lty = 1,bty = 'n',lwd = 2,seg.len = .2,x.intersp = 0.1,y.intersp = 0.5)
### Plot of VE over time
a3 <- ceiling(max(Nsusc)) 

VEmax <- max(VEls,trueVEls)

y3tck2 <- seq(0,a3,a3/5)/a3*VEmax
y3tck2lab <- sapply(round(seq(0,a3,a3/5)), function(x) toString(x))

# par(mar = c(5,4,4,5))

fpath <- paste0(bfolder,'Waning-Immunity-artefact/WIpresentation/WIgraphs')
setwd(fpath)

pdf('VE_2_virus.pdf',paper='USr',height = 8.5,width = 11) 
plot(times1[-c(1:100)],VEls[-c(1:100)],type = 'l',ylim = c(min(VEls,0),VEmax),ylab = 'VE',xlab = 'Day',
     main = '',lwd = 2)
lines(times1[-c(1:100)],trueVEls[-c(1:100)], col = 'red')

legend('left',c('Measured VE','True VE'),col = c('black','red'),
       lty = 1,bty = 'n',lwd = 2)
dev.off()

thetals <- 1 - cases1/(cases1 + cases2)
pdf('theta_time.pdf',paper='USr',height = 8.5,width = 11) 
plot(times1,thetals,type = 'l',ylab = expression(theta),xlab = 'Days',col = 'blue',lwd = 2,ylim = c(0,max(thetals)))
dev.off()
