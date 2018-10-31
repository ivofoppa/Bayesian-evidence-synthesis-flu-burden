###################################################################################################
###################################################################################################
Npop <- 1e+6

vacc <- .4
phi <- .6
csfrac <- 0.001 ### case sampling fraction
cnfrac <- csfrac*.1

lambda0 <- -log(0.85) ### 15% seasonal attack rate

xls <- seq(-3.5,3.5,length.out = 150)

lambdat0 <- dnorm(xls)
lambdat0 <- lambdat0 - min(lambdat0)
lambdat0 <- lambdat0/sum(lambdat0)*lambda0

lambdat1 <- lambdat0*(1-phi)


nsim <- 10000
VEpopls <- NULL

for (Npop in rev(seq(50,2000,50))){ ### pom: "population order of magnitude"
  
  VEls <- NULL
  sim <- 0
  while (sim < nsim){
    casetimels0 <- NULL
    N00 <- rbinom(1,Npop,vacc)
    N0 <- N00
    
    for (k in seq_along(lambdat0)){
      pinf <- 1 - exp(-lambdat0[k])
      ct <- rbinom(1,N0,pinf)
      casetimels0 <- c(casetimels0,ct)
      N0 <- N0 - ct
    }
    
    pt0_atrisk <- N0*length(xls) + sum(sapply(seq_along(casetimels0), function(t) casetimels0[t]*t))
    
    casetimels1 <- NULL
    N1 <- Npop - N00
    
    for (k in seq_along(lambdat1)){
      pinf <- 1 - exp(-lambdat1[k])
      ct <- rbinom(1,N1,pinf)
      casetimels1 <- c(casetimels1,ct)
      N1 <- N1 - ct
    }
    
    pt1_atrisk <- N1*length(xls) + sum(sapply(seq_along(casetimels1), function(t) casetimels1[t]*t))
    
    VEel <- 1 - sum(casetimels1)/pt1_atrisk/(sum(casetimels0)/pt0_atrisk)
    
    if (VEel!=-Inf){
      VEls <- c(VEls,VEel)
    }
    
    sim <- sim + 1
  }
  
  VEpopls <- c(VEpopls,mean(VEls[which(VEls!=Inf & VEls!=-Inf)]),ylab='Mean(VE estimate)',xlab = 'Source pop.')
  cat(paste0('Npop=',Npop,'\n'))
}

save(VEpopls, lambdat0, file = 'VEpopls.RData')

plot(tail(rev(seq(50,2000,50)),-1),VEpopls, type = 'l', ylim = c(.4,.65),ylab='Mean(VE estimate)',xlab = 'Source pop.')
lines(tail(rev(seq(50,2000,50)),-1),rep(.6,length(seq(50,2000,50)[-1])),col = 'red')


lambda0*(1-phi)/150

nsim <- 1000000
N1 <- rbinom(nsim,Npop,vacc)
N0 <- Npop - N1

phi <- .6 ### VE
lambda1 <- lambda0*(1-phi)

cases1 <- rpois(nsim,N1*lambda1*csfrac)
cases0 <- rpois(nsim,N0*lambda0*csfrac)

controls1 <- rpois(nsim,N1*cnfrac)
controls0 <- rpois(nsim,N0*cnfrac)

mean(1 - cases1/cases0/(controls1/controls0))

t1ls <- sapply(seq_along(cases1),function(k) (N1[k] - cases1[k]) + sum(rexp(cases1[k]),lambda1)
inft0ls <- rexp(cases0,lambda0)
inft1ls[1:10]
