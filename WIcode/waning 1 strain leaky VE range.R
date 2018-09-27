## test
library(deSolve)
bfolder <- "C:/Users/vor1/Dropbox/Misc work/Waning immunity/WI git project/"
setwd(paste0(bfolder,"WIwriteup"))

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
seas <- 365

dt <- .01
times <- seq(0, seas, by = dt)
###################################################################################################
###################################################################################################
v <- -log(1 - (vacc1 - vacc0))/seas
v <- 0
VErnge <- seq(0.2,.5,.1)

casearr <- VEarr <- trueVEarr <- infarr <- list()
for (k in seq_along(VErnge)){
  phi <- VErnge[k]
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
  
  rls <- casesv/casesnv
  vacc <- ((xv + yv + zv)/(xv + yv + zv +xnv + ynv + znv))[-1]
  vaccsusc <- (xv / (xv + xnv))[-1]
  
  casels <- casesnv + casesv
  infls <- yv + ynv
  VEls <- sapply(1:length(casesv), function(x) ifelse(casesv[x]<1e-5,NA,1-(rls[x]/(vacc[x]/(1-vacc[x])))))
  trueVEls <- sapply(1:length(casesv), function(x) ifelse(casesv[x]<1e-5,NA,1-(rls[x]/(vaccsusc[x]/(1-vaccsusc[x])))))
  # trueVEls <- sapply(VEls,function(x) phi)
  
  casearr[[k]] <- casels
  infarr[[k]] <- infls
  VEarr[[k]] <- VEls
  trueVEarr[[k]] <- trueVEls
}
###################################################################################################
###  Plots ########################################################################################
###  Trajectories are selected to start at the beginning of the epidemic   ########################
###################################################################################################
crit <- 10*dt
minsells <- maxsells <- NULL
for (k in seq_along(VErnge)){
  infls <- casearr[[k]]
  minsells <- c(minsells,min(which(infls > crit)))
  maxsells <- c(maxsells,max(which(infls > crit)))
}

# for (k in seq_along(VErnge)){
# 
#   selind <- unique(c(1:minsells[k],maxsells[k]:length(times)))
#   casearr[[k]][selind] <- NA
#   VEarr[[k]][selind] <- NA
# }
maxtime <- max(sapply(seq_along(VErnge), function(x) maxsells[x] - minsells[x]))

timesselList <- lapply(seq_along(VErnge), function(x) seq(minsells[x],minsells[x] + maxtime))
###################################################################################################
###  Saving workspace for use in Markdown document ################################################
###################################################################################################
filepath <- 'C:/Users/IFoppa/Documents/GitHub/Waning-Immunity-artefact/WIwriteup/WIplots/workspace.RData'
save.image(file = filepath)
# load(filepath)
###################################################################################################
###################################################################################################cols <- rainbow(length(VErnge))
### Epi curves
cols <- rainbow(length(VErnge))
setwd('C:/Users/IFoppa/Documents/GitHub/Waning-Immunity-artefact/WIwriteup/WIplots')

pdf('Epicurves.pdf',paper='USr',height = 8.5,width = 11) 
plot((0:maxtime)*dt,casearr[[1]][timesselList[[1]]],type = 'l', col = cols[1], ylab = 'Incidence', xlab = 'Day')

for (k in seq_along(VErnge)){
  lines((0:maxtime)*dt,casearr[[k]][timesselList[[1]]],col = cols[k])
}
legend('top',c('0.2','0.3','0.4','0.5'), lty = 1, col = cols,bty = 'n' ,title = 'VE')
dev.off()
###################################################################################################
### VE over time
pdf('VEtime.pdf',paper='USr',height = 8.5,width = 11) 
plot((0:maxtime)*dt, VEarr[[1]][timesselList[[1]]],type = 'l', col = cols[1], ylim = c(0,max(VErnge)*1.1), ylab = 'VE est.', xlab = 'Day')

for (k in seq_along(VErnge)){
  lines((0:maxtime)*dt, VEarr[[k]][timesselList[[k]]],col = cols[k])
}
legend('top',rev(c('0.2','0.3','0.4','0.5')), lty = 1, col = rev(cols),bty = 'n' ,title = 'VE')
dev.off()
####################################################################################################
####################################################################################################
### Abs. VE bias over time
pdf('VEbias_abs.pdf',paper='USr',height = 8.5,width = 11) 

miny <- min((VEarr[[1]][timesselList[[1]]] - VErnge[1]),na.rm = T)
maxy <- max((VEarr[[1]][timesselList[[1]]] - VErnge[1]),na.rm = T)

plot((0:maxtime)*dt, (VEarr[[1]][timesselList[[1]]] - VErnge[1]),type = 'l', col = cols[1], ylab = 'Abs. Bias', xlab = 'Day',yaxt = 'n',ylim = c(-.22,0))
axis(2,c(-.2,-.15,-.10,-.05,0),labels = c(-20,-15,-10,-5,0))

for (k in seq_along(VErnge)){
  lines((0:maxtime)*dt, (VEarr[[k]][timesselList[[k]]] - VErnge[k]),col = cols[k])
}
legend('left',rev(c('0.2','0.3','0.4','0.5')), lty = 1, col = rev(cols),bty = 'n' ,title = 'VE')
dev.off()
####################################################################################################
####################################################################################################
### Rel. VE bias over time
miny <- min((VEarr[[1]][timesselList[[1]]] - VErnge[1])/VErnge[1],na.rm = T)

pdf('VEbias_rel.pdf',paper='USr',height = 8.5,width = 11) 

plot((0:maxtime)*dt, (VEarr[[1]][timesselList[[1]]] - VErnge[1])/VErnge[1],type = 'l', col = cols[1], ylab = 'Rel. Bias', xlab = 'Day',yaxt = 'n')
axis(2,c(-.75,-.5,-.25,0),labels = c(-75,-50,-25,0))

for (k in seq_along(VErnge)){
  lines((0:maxtime)*dt, (VEarr[[k]][timesselList[[k]]] - VErnge[k])/VErnge[k],col = cols[k])
}
legend('left',rev(c('0.2','0.3','0.4','0.5')), lty = 1, col = rev(cols),bty = 'n' ,title = 'VE')
dev.off()
####################################################################################################
####################################################################################################
### Attack rate
pdf('AR.pdf',paper='USr',height = 8.5,width = 11) 
plot((0:maxtime)*dt, cumsum(casearr[[1]][timesselList[[1]]])/Ntot,type = 'l', col = cols[1], ylim = c(0,1), ylab = 'Attack rate', xlab = 'Day')

for (k in seq_along(VErnge)){
  lines((0:maxtime)*dt, cumsum(casearr[[k]][timesselList[[k]]])/Ntot,col = cols[k])
}
legend('top',c('0.2','0.3','0.4','0.5'), lty = 1, col = cols,bty = 'n' ,title = 'VE')
dev.off()
###################################################################################################