library(deSolve)
bfolder <- "C:/Users/vor1/Dropbox/Misc work/Waning immunity/WI project/"
setwd(paste0(bfolder,"WIgraphs"))

## asign values to parameters
Ntot <- 2000000
vacc0 <- 0.38
# vacc1 <- 0.47
vacc1 <- vacc0
##################################################################################################
##  Initial values ################################################################################
###################################################################################################
###################################################################################################
Nv0 <- round(Ntot*vacc)
Nnv0 <- Ntot - Nv0

gamma <- 1/4

beta <- 2*gamma/(Ntot)
epsilon <- 0.2 # proportion preexisting immunity
### Initial conditions
yv0 <- yv <- 0; 

ynv0 <- ynv <- 1; 

z0 <- 0
###################################################################################################
###################################################################################################
seas <- 150

dt <- 0.01
times <- seq(0, seas, by = dt)

nsims <- length(times) - 1
###################################################################################################
###################################################################################################
casearr <- array(0,dim = c(0,nsims))
VEarr <- array(0,dim = c(0,nsims))
trueVEarr <- array(0,dim = c(0,nsims))

v <- (vacc1 - vacc0)/seas*dt

VErnge <- c(0.2,.5,.7)
for (phi in VErnge){
  xv0 <- Nv0 * (1 - epsilon)*(1-phi) - yv0
  xnv0 <- Nnv0 * (1 - epsilon) - ynv0
  ### define parameter vector for ode function
  parameters <- c(gamma=gamma,beta=beta,phi=phi)
  ### Define model
  KKmod <- function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      # rate of change
      dxv <- - beta*(1-phi)*(1-v)*xv*(yv + ynv) 
      dxnv <- - beta*xnv*(1-v)*(yv + ynv)
      
      dyv <- beta*(1-phi)*xv*(yv + ynv) - yv*gamma
      dynv <- beta*xnv*(yv + ynv) - ynv*gamma

      dz <- (ynv + yv)*gamma
      
      # return the rate of change
      list(c(dxv,dxnv,dynv,dyv,dz))
    })
  }
  ### define vector with initial conditions for ode function
  state <- c(xv=xv0,xnv=xnv0,ynv=ynv0,yv=yv0,z=z0)
  

  out <- data.frame(ode(y = state, times = times, func = KKmod, parms = parameters))
  
  xv <- out$xv
  xnv <- out$xnv
  yv <- out$yv
  ynv <- out$ynv
  
  ###################################################################################################
  ###################################################################################################
  casesv <- -diff(xv)
  casesnv <- -diff(xnv)
  
  Nsusc <- (xv + xnv)[-1]
  
  infprsv <- (yv + ynv)[-1]*beta*(1-phi)
  infprsnv <- (yv + ynv)[-1]*beta
  
  rls <- casesv/casesnv
  
  casels <- casesnv + casesv
  VEls <- 1-(rls/(vacc/(1-vacc)))
  trueVEls <- 1-infprsv/infprsnv
  
  casearr <- rbind(casearr,casels,deparse.level = 0)
  VEarr <- rbind(VEarr,VEls,deparse.level = 0)
  trueVEarr <- rbind(trueVEarr,trueVEls,deparse.level = 0)
}
###################################################################################################
###################################################################################################