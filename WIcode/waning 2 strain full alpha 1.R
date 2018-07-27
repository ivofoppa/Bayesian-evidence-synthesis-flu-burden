### Investigating the presence of waning when all-or none effect with against specific influenza strains
library(deSolve)
## asign values to parameters
Ntot <- 20000000
vacc <- 0.3
# phi <- 1
Nv0 <- round(Ntot*vacc)
Nnv0 <- Ntot - Nv0
alpha <- 1 #VE

Nv10 <- Nv0*alpha
Nv20 <- Nv0-Nv10
gamma <- 0.333

beta1 <- 1.25*gamma/(Ntot)
beta2 <- 1.75*gamma/(Ntot)
### Initial conditions
yv110 <- yv11 <- 2/3/5; 
yv210 <- yv21 <- 2/3/5; 
yv220 <- yv22 <- 2/3/5; 

ynv10 <- ynv1 <- 1; 
ynv20 <- ynv2 <- 1; 

xv10 <- Nv10 - yv110
xv20 <- Nv20 - yv210 - yv210
xnv0 <- Nnv0 - ynv10 - ynv20
z0 <- 0
parameters <- c(gamma=gamma,beta1=beta1, beta2=beta2)
### Define model
KKmod <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dxv1 <- - beta1*xv1*(yv11 + yv21 + ynv1) 
    dxv2 <- - beta1*xv2*(yv11 + yv21 + ynv1) - beta2*xv2*(yv22 + ynv2) 
    dxnv <- - beta1*xnv*(yv11 + yv21 + ynv1) - beta2*xnv*(yv22 + ynv2)
    
    dyv11 <- beta1*xv1*(yv11 + yv21 + ynv1) - yv11*gamma
    
    dyv21 <- beta1*xv2*(yv11 + yv21 + ynv1) - yv21*gamma
    dyv22 <- beta2*xv2*(yv22 + ynv2) - yv22*gamma
    
    dynv1 <- beta1*xnv*(yv11 + yv21 + ynv1) - ynv1*gamma
    dynv2 <- beta2*xnv*(yv22 + ynv2) - ynv2*gamma
    
    dz <- (ynv1 + ynv2 + yv11 + yv21 + yv22)*gamma
    
    # return the rate of change
    list(c(dxv1,dxv2,dxnv,dynv1,dynv2,dyv11,dyv21,dyv22, dz))
  })
}
### define vector with initial conditions for ode function
state <- c(xv1=xv10,xv2=xv20,xnv=xnv0,ynv1=ynv10,ynv2=ynv20,yv11=yv110,yv21=yv210,yv22=yv220,z=z0)

times <- seq(0, 300, by = 0.01)

out <- data.frame(ode(y = state, times = times, func = KKmod, parms = parameters))

xv1 <- out$xv1
xv2 <- out$xv2
xnv <- out$xnv
yv11 <- out$yv11
yv21 <- out$yv21
yv22 <- out$yv22
ynv1 <- out$ynv1
ynv2 <- out$ynv2
###################################################################################################
###################################################################################################

# dev.off()

