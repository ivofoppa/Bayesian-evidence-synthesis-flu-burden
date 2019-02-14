library(deSolve)
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Waning-Immunity-artefact/' ## Define root path (where repository)
## asign values to parameters
Ntot <- 2000000
vacc0 <- 0.38
vacc1 <- 0.47
vacc0 <- vacc1
##################################################################################################
##  Initial values ################################################################################
###################################################################################################
###################################################################################################
Nv0 <- round(Ntot*vacc0)
Nnv0 <- Ntot - Nv0

gamma <- 1/4
r0 <- 1.6
beta <- r0*gamma/(Ntot)
epsilon <- 0 # proportion preexisting immunity
###################################################################################################
seas <- 500

dt <- .01
times <- seq(0, seas, by = dt)
###################################################################################################
###################################################################################################
v <- -log(1 - (vacc1 - vacc0))/seas
v <- 0
VErnge <- seq(0.2,.6,.1)

vk <- 4 ## VE index: VE= 0.5

phi <- VErnge[vk]
###################################################################################################
### define parameter vector for ode function
parameters <- c(gamma=gamma,beta=beta,phi=phi,v=v)
### Define model
KKmod <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dxv <- - beta*(1-phi)*xv*(yv + ynv) + xnv*v
    dxnv <- - beta*xnv*(yv + ynv) - v*xnv
    
    dcasesv <- beta*(1-phi)*xv*(yv + ynv)
    dcasesnv <- beta*xnv*(yv + ynv)
    
    dyv <- beta*(1-phi)*xv*(yv + ynv) - yv*gamma
    dynv <- beta*xnv*(yv + ynv) - ynv*gamma
    
    dzv <- yv*gamma
    dznv <- ynv*gamma
    
    # return the rate of change
    list(c(dxv,dxnv,dynv,dyv,dzv,dznv,dcasesv,dcasesnv))
  })
}
### Initial conditions
y0 <- 1
denom <- vacc0*(1-phi) + (1-vacc0)
yv0 <- vacc0*(1-phi)/denom*y0; 
ynv0 <- y0 - yv0; 

xv0 <- Nv0 * (1 - epsilon) - yv0
xnv0 <- Nnv0 * (1 - epsilon) - ynv0

zv0 <- znv0 <- 0
###################################################################################################
state <- c(xv=xv0,xnv=xnv0,ynv=ynv0,yv=yv0,zv=0,znv=0,casesv=0,casesnv=0)
out <- data.frame(ode(y = state, times = times, func = KKmod, parms = parameters))

xv <- out$xv
xnv <- out$xnv
yv <- out$yv
ynv <- out$ynv
zv <- out$zv
znv <- out$znv
casesv <- out$casesv
casesnv <- out$casesnv
###################################################################################################
###################################################################################################
casesv <- diff(casesv)
casesnv <- diff(casesnv)

case_v <- sapply(1:seas, function(t) round(sum(casesv[which(times[-1] > (t - 1) & times[-1] <= t)])))
case_nv <- sapply(1:seas, function(t) round(sum(casesnv[which(times[-1] > (t - 1) & times[-1] <= t)])))
data <- data.frame('vacc'=case_v,'unvacc'=case_nv)

setwd(bfolder)
outfname <- 'simuldata.csv'
write.csv(data,outfname,row.names = FALSE)
###################################################################################################
###################################################################################################
###################################################################################################