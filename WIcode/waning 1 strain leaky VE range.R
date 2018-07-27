library(deSolve)
bfolder <- "C:/Users/vor1/Dropbox/Misc work/Waning immunity/WI project/"
setwd(paste0(bfolder,"WIgraphs"))

## asign values to parameters
Ntot <- 2000000
vacc <- 0.3
##################################################################################################
##  Initial values ################################################################################
###################################################################################################
###################################################################################################
Nv0 <- round(Ntot*vacc)
Nnv0 <- Ntot - Nv0

gamma <- 0.333

beta <- 1.8*gamma/(Ntot)
### Initial conditions
yv0 <- yv <- 0; 

ynv0 <- ynv <- 1; 

xv0 <- Nv0 - yv0
xnv0 <- Nnv0 - ynv0
z0 <- 0
###################################################################################################
###################################################################################################
seas <- 150
times <- seq(0, seas, by = 0.01)

nsims <- length(times) - 1
###################################################################################################
###################################################################################################
casearr <- array(0,dim = c(0,nsims))
VEarr <- array(0,dim = c(0,nsims))
trueVEarr <- array(0,dim = c(0,nsims))

for (phi in seq(0.1,.9,.1)){
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
  

  out <- data.frame(ode(y = state, times = times, func = KKmod, parms = parameters))
  
  xv <- out$xv
  xnv <- out$xnv
  
  ###################################################################################################
  ###################################################################################################
  dxv <- -diff(xv)
  dxnv <- -diff(xnv)

  casesv <- dxv
  casesnv <- dxnv
  
  Nsusc <- (xv + xnv)[-1]
  
  infprsv <- (out$yv + out$ynv)[-1]*beta*(1-phi)
  infprsnv <- (out$yv + out$ynv)[-1]*beta
  
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