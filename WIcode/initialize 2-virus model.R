library(binhf)
library(survival)
Ntot <- 2000000

pSI <- function(r01,I1,r02,I2){ ## infection probability
  beta1 <- r01 * delta
  beta2 <- r02 * delta
  if((beta1 * I1 + beta2 * I2)==0){
    return (c(0,0,1))
  } else {
    lambda <- (beta1 * I1 + beta2 * I2) / Ntot
    pinf <- as.numeric(1-exp(-lambda))
    p1 <- beta1*I1/(beta1*I1 + beta2*I2)
    p2 <- 1 - p1
    return (c(pinf*p1,pinf*p2,1 - pinf))
  }
}

vdurls <- 1:vdur
if (vacc1 !=0){
  raw_vls <- sapply(vdurls,function(x) dnorm(x,vdur/2,vdur/5)) ## Normally distributed vacc uptake
  s1 <- -log(vacc1)
  v0 <- raw_vls/sum(raw_vls) * s1 ## scaled right
  v <- c(v0,rep(0,seas + prevdur - vdur))
} else {
  v <- c(rep(0,seas + prevdur))
}
###################################################################################################
prem <- 1 - exp(-delta) ## Daily removal probability
###################################################################################################
vaccdelim <- unique(c(seq(1,seas + prevdur,vaccint),seas + prevdur))
nvaccat <- length(vaccdelim) - 1 ## number of time-since-vacc categories
###################################################################################################
Snvinit <- rmultinom(1,Ntot,c(vs1,vs2,vs0,vs12))
S1nvinit <- Snvinit[1] ## only susceptible to virus 1 if vacc.
S2nvinit <- Snvinit[2] ## only susceptible to virus 2 if vacc.
S0nvinit <- Snvinit[3] ## susceptible to neither virus if vacc.
S12nvinit <- Snvinit[4] ## remain susceptible to either virus if vacc.

S1vinit <- S2vinit <- S0vinit <- S12vinit <-  rep(0,seas + prevdur) ## Susceptible to only virus 1, 2 neither if vacc.

I1nvinit <- I2nvinit <- 0
I1vinit <- I2vinit <- rep(0,seas + prevdur)

Rnvinit <- 0
Rvinit <- rep(0,seas + prevdur)

inits <- list(S1vinit=S1vinit,S2vinit=S2vinit,S0vinit=S0vinit,S12vinit=S12vinit,
              S1nvinit=S1nvinit,S2nvinit=S2nvinit,S0nvinit=S0nvinit,S12nvinit=S12nvinit,
              I1vinit=I1vinit,I2vinit=I2vinit,
              I1nvinit=I1nvinit,I2nvinit=I2nvinit,
              Rvinit=Rvinit,Rnvinit=Rnvinit)
###################################################################################################
###################################################################################################
###################################################################################################
S1nv <- as.integer(inits$S1nvinit)
S2nv <- as.integer(inits$S2nvinit)
S0nv <- as.integer(inits$S0nvinit)
S12nv <- as.integer(inits$S12nvinit)

S1v <- as.integer(inits$S1vinit)
S2v <- as.integer(inits$S2vinit)
S0v <- as.integer(inits$S0vinit) 
S12v <- as.integer(inits$S12vinit)

I1nv <- as.integer(inits$I1nvinit)
I2nv <- as.integer(inits$I2nvinit)

I1v <- as.integer(inits$I1vinit)
I2v <- as.integer(inits$I2vinit)

Rnv <- as.integer(inits$Rnvinit)
Rv <- as.integer(inits$Rvinit)

time <- 1

while (time <= prevdur){
  ###################################################################################################
  ###################################################################################################
  ### Vaccination
  pvacc <- 1 - exp(-v[time]) ## vaccination uptake at this time
  
  S1v <- shift(S1v,1)
  s1vacc <- rbinom(1,S1nv,pvacc)
  S1v[1] <- s1vacc; S1nv <- S1nv - s1vacc
  
  S2v <- shift(S2v,1)
  s2vacc <- rbinom(1,S2nv,pvacc)
  S2v[1] <- s2vacc; S2nv <- S2nv - s2vacc
  
  S0v <- shift(S0v,1)
  s0vacc <- rbinom(1,S0nv,pvacc)
  S0v[1] <- s0vacc; S0nv <- S0nv - s0vacc
  
  S12v <- shift(S12v,1)
  s12vacc <- rbinom(1,S12nv,pvacc)
  S12v[1] <- s12vacc; S12nv <- S12nv - s12vacc
  
  I1v <- shift(I1v,1) ### Infectious not getting vaccinated, by assumption
  I2v <- shift(I2v,1) ### note that number represent virus, not type
  # 
  Rv <- shift(Rv,1)
  rvacc <- rbinom(1,Rnv,pvacc)
  Rv[1] <- rvacc; Rnv <- Rnv - rvacc
  
  time <- time + 1
}

## Infections in unvaccinated
pinf11nv <-  S1nv/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v)) ### probability that infections with virus 1 are in type 1 individuals
pinf21nv <-  S2nv/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v))  ### probability that infections with virus 1 are in type 2 individuals
pinf01nv <-  S0nv/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v))  ### probability that infections with virus 1 are in type 0 individuals
pinf121nv <-  S12nv/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v))  ### probability that infections with virus 1 are in type 12 individuals
## virus 2
pinf12nv <-  S1nv/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v)) ### probability that infections with virus 1 are in type 1 individuals
pinf22nv <-  S2nv/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v))  ### probability that infections with virus 1 are in type 2 individuals
pinf02nv <-  S0nv/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v))  ### probability that infections with virus 1 are in type 0 individuals
pinf122nv <-  S12nv/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v))  ### probability that infections with virus 1 are in type 12 individuals
## Infections in vaccinated
pinf11v <-  sum(S1v)/(S1nv + S2nv + S0nv + S12nv + sum(S1v) + sum(S12v))
pinf121v <-  1 - pinf11nv - pinf21nv - pinf01nv - pinf121nv - pinf11v  ### probability that infections with virus 1 are in type 1 individuals
## virus 2
pinf22v <-  sum(S2v)/(S1nv + S2nv + S0nv + S12nv + sum(S2v) + sum(S12v))
pinf122v <-  1 - pinf12nv - pinf22nv - pinf02nv - pinf122nv - pinf22v  ### probability that infections with virus 1 are in type 1 individuals

## infection numbers by virus and type
## virus 1
inf1 <- rmultinom(1,inf1num,c(pinf11nv,pinf21nv,pinf01nv,pinf121nv,pinf11v,pinf121v))
inf11nv <- inf1[1]
inf21nv <- inf1[2]
inf01nv <- inf1[3]
inf121nv <- inf1[4]
inf11vnum <- inf1[5]
inf11v <- c(rmultinom(1,inf11vnum,rep(1/prevdur,prevdur)),rep(0,seas))
inf121vnum <- inf1[6]
inf121v <- c(rmultinom(1,inf121vnum,rep(1/prevdur,prevdur)),rep(0,seas))

## virus 2
inf2 <- rmultinom(1,inf2num,c(pinf12nv,pinf22nv,pinf02nv,pinf122nv,pinf22v,pinf122v))
inf12nv <- inf2[1]
inf22nv <- inf2[2]
inf02nv <- inf2[3]
inf122nv <- inf2[4]
inf22vnum <- inf2[5]
inf22v <- c(rmultinom(1,inf22vnum,rep(1/prevdur,prevdur)),rep(0,seas))
inf122vnum <- inf2[6]
inf122v <- c(rmultinom(1,inf122vnum,rep(1/prevdur,prevdur)),rep(0,seas))

#########################################################################################
#########################################################################################
## updating states: Infections with virus 1
I1nv <- I1nv + inf11nv + inf21nv + inf01nv + inf121nv
I1v <- I1v + inf11v + inf121v
## updating states: Infections with virus 2
I2nv <- I2nv + inf12nv + inf22nv + inf02nv + inf122nv
I2v <- I2v + inf22v + inf122v

## Updating # susceptibles
## unvaccinated
S1nv <- S1nv - inf11nv - inf12nv
S2nv <- S2nv - inf21nv - inf22nv
S0nv <- S0nv - inf01nv - inf02nv
S12nv <- S12nv - inf121nv - inf122nv
## vaccinated
S1v <- S1v - inf11v
S2v <- S2v - inf22v
S12v <- S12v - inf121v - inf122v

I1 <- sum(c(I1nv,I1v))
I2 <- sum(c(I2nv,I2v))

I1ls <- I2ls <- kls <- NULL
