library(deSolve)
bfolder <- "C:/Users/vor1/Dropbox/Misc work/Waning immunity/WIpresentations/"
setwd(paste0(bfolder,"WIgraphs"))

## asign values to parameters
Ntot <- 20000000
vacc <- 0.5
phi <- .5
Nv0 <- round(Ntot*vacc)
Nnv0 <- Ntot - Nv0

seas <- 300

gamma <- 0.333
###################################################################################################
###  Dynamic model ################################################################################
###################################################################################################
beta <- 1.8*gamma/(Ntot)
### Initial conditions
yv0 <- yv <- 0; 

ynv0 <- ynv <- 1; 

xv0 <- Nv0 - yv0
xnv0 <- Nnv0 - ynv0
z0 <- 0
### define parameter vector for ode function
parameters <- c(gamma=gamma,beta=beta,phi=phi)
### Define model
KKmod <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dxv <- - beta*(1-phi)*xv*(yv + ynv) 
    dxnv <- - beta*xnv*(yv + ynv)
    
    dyv <- beta*(1-phi)*xv*(yv + ynv) - yv*gamma
    dynv <- beta*xnv*(yv + ynv) - ynv*gamma
    
    dyvnew <- beta*(1-phi)*xv*(yv + ynv)
    dynvnew <- beta*xnv*(yv + ynv)
    
    dz <- (ynv + yv)*gamma
    
    # return the rate of change
    list(c(dxv,dxnv,dynv,dyv,
           dynvnew,dyvnew, dz))
  })
}
### define vector with initial conditions for ode function
state <- c(xv=xv0,xnv=xnv0,ynv=ynv0,yv=yv0,
           ynvnew=0,yvnew=0,z=z0)

times <- seq(0, seas, by = 0.01)

out <- data.frame(ode(y = state, times = times, func = KKmod, parms = parameters))

xv <- out$xv
xnv <- out$xnv
yvnew <- out$yvnew
ynvnew <- out$ynvnew
###################################################################################################
###################################################################################################
dxv <- -diff(xv)
dxnv <- -diff(xnv)

dynv <- diff(ynvnew)
dyv <- diff(yvnew)

selind <- which((dxv1 + dxv2 + dxnv)>0)

casesv <- dyv; casesv <- casesv[selind]
casesnv <- dynv; casesnv <- casesnv[selind]

newcases <- casesv + casesnv;

Nsusc <- (xv + xnv)[-1]

infprsv <- (dyv + dynv)*beta*(1-phi)
infprsnv <- (dyv + dynv)*beta

rls <- casesv/casesnv

VEls <- 1-( rls/(vacc/(1-vacc)))
trueVEls <- 1-infprsv/infprsnv

### Plot of VE over time
a3 <- ceiling(max(Nsusc)) 
b3 <- nchar(format(a3,scientific = F))
trueVEls <- 1-infprsv/infprsnv

y3uplim0 <- ceiling(a3/10^(b3-2))
y3uplim <- ifelse(y3uplim0%%2==0,y3uplim0,y3uplim0+1)/100

y3tck2 <- seq(0,y3uplim,y3uplim/4)*10^b3/a3*VEmax
y3tck2lab <- sapply(round(seq(0,y3uplim,y3uplim/4)*10^b3), function(x) toString(x))

y3toax <- y3uplim0*10^(b3-1)

par(mar = c(5,4,4,5))

VEmax <- max(VEls)
plot(times1,VEls,type = 'l',ylim = c(min(VEls,0),VEmax),ylab = 'VE',xlab = 'Day',
     main = '',lwd = 2)
lines(times1,trueVEls, col = 'red')
lines(times1,Nsusc/a3*VEmax, col = 'blue',lwd = 2)
axis(side = 4,at= y3tck2,labels = y3tck2lab)
mtext('Susceptibles', 4, line = 3)

legend('bottomleft',c('VE','True VE','Pop. susc.'),col = c('black','red','blue'),
       lty = 1,bty = 'n',lwd = 2,seg.len = 1)

