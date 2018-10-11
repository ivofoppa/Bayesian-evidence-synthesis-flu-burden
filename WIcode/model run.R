studydata0 <- array(0,dim = c(seas,1 + 2 * nvaccat + 2))
studydata10 <- array(0,dim = c(seas,1 + 2 * nvaccat + 2))
studydata20 <- array(0,dim = c(seas,1 + 2 * nvaccat + 2))

while (time <= seas + prevdur){
  pinf12 <- pSI(r10,I1,r20,I2) ## for fully susceptible
  pinf1 <- pSI(r10,I1,r20,0) ## for type 1 vacc.
  pinf2 <- pSI(r10,0,r20,I2) ## for type 2 vacc.
  ## unvaccinated
  inf1nvls <- rmultinom(1, S1nv, pinf12)
  inf2nvls <- rmultinom(1, S2nv, pinf12)
  inf0nvls <- rmultinom(1, S0nv, pinf12)
  inf12nvls <- rmultinom(1, S12nv, pinf12)
  ## vaccinated
  inf1varr <- sapply(S1v, function(s1v) rmultinom(1, s1v, pinf1))
  inf2varr <- sapply(S2v, function(s2v) rmultinom(1, s2v, pinf2))
  inf12varr <- sapply(S12v, function(s12v) rmultinom(1, s12v, pinf12))
  
  ## unvaccinated
  ### virus 1
  inf11nv <- inf1nvls[1]
  inf21nv <- inf2nvls[1]
  inf01nv <- inf0nvls[1]
  inf121nv <- inf12nvls[1]
  ### virus 2
  inf12nv <- inf1nvls[2]
  inf22nv <- inf2nvls[2]
  inf02nv <- inf0nvls[2]
  inf122nv <- inf12nvls[2]
  ### updating susceptibles
  S1nv <- inf1nvls[3]
  S2nv <- inf2nvls[3]
  S0nv <- inf0nvls[3]
  S12nv <- inf12nvls[3]
  
  ## vaccinated
  ### virus 1
  inf11v <- inf1varr[1,]
  inf11v2 <- sapply(seq_along(inf11v), function(k) ifelse(k <= (time - prevdur), 0, inf11v[k]))
  
  inf121v <- inf12varr[1,]
  inf121v2 <- sapply(seq_along(inf121v), function(k) ifelse(k <= (time - prevdur), 0, inf121v[k]))
  
  ### virus 2
  inf22v <- inf2varr[2,]
  inf22v2 <- sapply(seq_along(inf22v), function(k) ifelse(k <= (time - prevdur), 0, inf22v[k]))
  
  inf122v <- inf12varr[2,]
  inf122v2 <- sapply(seq_along(inf122v), function(k) ifelse(k <= (time - prevdur), 0, inf122v[k]))
  
  ### updating susceptibles
  S1v <- inf1varr[3,]
  S2v <- inf2varr[3,]
  S12v <- inf12varr[3,]

  ## for plotting purposes ...
  I1ls <- c(I1ls,sum(inf11v) + sum(inf121v) + inf11nv + inf21nv + inf01nv + inf121nv)  
  ## virus 2
  I2ls <- c(I2ls,sum(inf22v) + sum(inf122v) + inf12nv + inf22nv + inf02nv + inf122nv)  
  ## list of "k factors'
  kfact <- (sum(S1v + S12v)*r10*delta*I1 + sum(S2v + S12v)*r20*delta*I2)*(S1nv + S2nv + S0nv + S12nv)/
    (sum(S1v + S2v + S12v))/
    ((S1nv + S2nv + S0nv + S12nv)*r10*delta*I1 + (S1nv + S2nv + S0nv + S12nv)*r20*delta*I2)
  kls <- c(kls,kfact)
  ## updating other states: Infections with virus 1
  I1nv <- I1nv + inf11nv + inf21nv + inf01nv + inf121nv
  I1v <- I1v + inf11v + inf121v
  ## updating other states: Infections with virus 2
  I2nv <- I2nv + inf12nv + inf22nv + inf02nv + inf122nv
  I2v <- I2v + inf22v + inf122v
  ## updating states: Infections with either virus
  ## numbers removed from virus 1 infectious
  rem1nv <- rbinom(1,I1nv,prem)
  rem1v <- sapply(I1v, function(iv) rbinom(1,iv, prem))
  ## numbers removed from virus 2 infectious
  rem2nv <- rbinom(1,I2nv,prem)
  rem2v <- sapply(I2v, function(iv) rbinom(1,iv, prem))
  ## updating infections with removals
  I1nv <- I1nv - rem1nv
  I1v <- I1v - rem1v
  I2nv <- I2nv- rem2nv
  I2v <- I2v - rem2v
  ## numbers removed from virus 2 infectious
  Rnv <- Rnv + rem1nv + rem2nv
  Rv <- Rv + rem1v + rem2v
  
  I1 <- sum(c(I1nv,I1v))
  I2 <- sum(c(I2nv,I2v))
  ### Interrupt evaluation if no more transmission
  if (I1==0 & I2==0 | is.na(I1) | is.na(I2)){
    break()
  }
  ###################################################################################################
  ### "Study" is conducted ##########################################################################
  ###################################################################################################
  infnv <- inf11nv + inf21nv + inf01nv + inf121nv + inf12nv + inf22nv + inf02nv + inf122nv
  infv <- inf11v + inf22v + inf121v + inf122v
  infv2 <- inf11v2 + inf22v2 + inf121v2 + inf122v2
  
  casesls <- c(infnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(infv[vaccdelim[x] : vaccdelim[x + 1]])))
  casesls2 <- c(infnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(infv2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  controlsvsum <- S1v + S2v + S0v + S12v + Rv
  controlsvsum1 <- S1v + S2v + S0v + S12v ### "True" comparison: only susceptibles
  controlsvsum2 <- sapply(seq_along(controlsvsum), function(k) ifelse(k <= (time - prevdur), 0, controlsvsum[k]))
  
  controlsnv <- S1nv + S2nv + S0nv + S12nv + Rnv
  controlsnv1 <- S1nv + S2nv + S0nv + S12nv
  
  controlsls <- c(controlsnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum[vaccdelim[x] : vaccdelim[x + 1]])))
  controlsls1 <- c(controlsnv1,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum1[vaccdelim[x] : vaccdelim[x + 1]])))
  controlsls2 <- c(controlsnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  studydata0[time - prevdur,] <- c(time - prevdur,casesls,controlsls)
  studydata10[time - prevdur,] <- c(time - prevdur,casesls,controlsls1)
  studydata20[time - prevdur,] <- c(time - prevdur,casesls2,controlsls2)
  ###################################################################################################
  ### Vaccination: proceeding time since vacc. and adding new vaccinees
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
  ###################################################################################################
  ### Stop when no more transmission 
  time <- time + 1
}

