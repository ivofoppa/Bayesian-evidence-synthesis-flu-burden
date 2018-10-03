source('initialize 2-virus model.R') ## run initialization code

## Running code to near break point
while (time <= 50 + prevdur){
  ## virus 1
  p1inf <- pSE(r10,I1)
  inf11nv <- rbinom(1, S1nv, p1inf)
  inf21nv <- rbinom(1, S2nv, p1inf)
  inf01nv <- rbinom(1, S0nv, p1inf)
  inf121nv <- rbinom(1, S12nv, p1inf)
  
  inf11v <- sapply(S1v, function(sv) rbinom(1, sv, p1inf))
  inf11v2 <- sapply(seq_along(inf11v), function(k) ifelse(k <= (time - prevdur), 0, inf11v[k]))
  inf121v <- sapply(S12v, function(sv) rbinom(1, sv, p1inf))
  inf121v2 <- sapply(seq_along(inf121v), function(k) ifelse(k <= (time - prevdur), 0, inf11v[k]))
  
  rem1nv <- rbinom(1,I1nv,prem)
  rem1v <- sapply(I1v, function(iv) rbinom(1,iv, prem))
  ## Updating # susceptibles, because used again for virus 2
  ## unvaccinated
  S1nv <- S1nv - inf11nv
  S2nv <- S2nv - inf21nv
  S0nv <- S0nv - inf01nv
  S12nv <- S12nv - inf121nv
  ## vaccinated
  S1v <- S1v - inf11v
  S12v <- S12v - inf121v
  
  ## for plotting purposes ...
  I1ls <- c(I1ls,sum(inf11v) + sum(inf121v) + inf11nv + inf21nv + inf01nv + inf121nv)  
  ## virus 2
  p2inf <- pSE(r20,I2)
  inf12nv <- rbinom(1, S1nv, p2inf)
  inf22nv <- rbinom(1, S2nv, p2inf)
  inf02nv <- rbinom(1, S0nv, p2inf)
  inf122nv <- rbinom(1, S12nv, p2inf)
  
  inf22v <- sapply(S2v, function(sv) rbinom(1, sv, p2inf))
  inf22v2 <- sapply(seq_along(inf22v), function(k) ifelse(k <= (time - prevdur), 0, inf22v[k]))
  inf122v <- sapply(S12v, function(sv) rbinom(1, sv, p1inf))
  inf122v2 <- sapply(seq_along(inf122v), function(k) ifelse(k <= (time - prevdur), 0, inf11v[k]))
  
  rem2nv <- rbinom(1,I2nv,prem)
  rem2v <- sapply(I2v, function(iv) rbinom(1,iv, prem))
  ## Updating # susceptibles
  ## unvaccinated
  S1nv <- S1nv - inf12nv
  S2nv <- S2nv - inf22nv
  S0nv <- S0nv - inf02nv
  S12nv <- S12nv - inf122nv
  ## vaccinated
  S2v <- S2v - inf22v
  S12v <- S12v - inf122v
  
  ## for plotting purposes ...
  I2ls <- c(I2ls,sum(inf22v) + sum(inf122v) + inf12nv + inf22nv + inf02nv + inf122nv)  
  #########################################################################################
  ## updating other states: Infections with virus 1
  I1nv <- I1nv + inf11nv + inf21nv + inf01nv + inf121nv
  I1v <- I1v + inf11v + inf121v
  ## updating other states: Infections with virus 2
  I2nv <- I2nv + inf12nv + inf22nv + inf02nv + inf122nv
  I2v <- I2v + inf22v + inf122v
  ## updating states: Infections with either virus
  Rnv <- Rnv + rem1nv + rem2nv
  Rv <- Rv + rem1v + rem2v
  
  I1 <- sum(c(I1nv,I1v))
  I2 <- sum(c(I2nv,I2v))
  ###################################################################################################
  ### "Study" is conducted ##########################################################################
  ###################################################################################################
  infnv <- inf11nv + inf21nv + inf01nv + inf121nv + inf12nv + inf22nv + inf02nv + inf122nv
  infv <- inf11v + inf22v + inf121v + inf122v
  infv2 <- inf11v2 + inf22v2 + inf121v2 + inf122v2
  
  casesls <- c(infnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(infv[vaccdelim[x] : vaccdelim[x + 1]])))
  casesls2 <- c(infnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(infv2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  controlsvsum <- S1v + S2v + S0v + S12v + Rv
  controlsvsum2 <- sapply(seq_along(controlsvsum), function(k) ifelse(k <= (time - prevdur), 0, controlsvsum[k]))
  
  controlsnv <- S1nv + S2nv + S0nv + S12nv + Rnv
  
  controlsls <- c(controlsnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum[vaccdelim[x] : vaccdelim[x + 1]])))
  controlsls2 <- c(controlsnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  studydata[time - prevdur,] <- c(time - prevdur,casesls,controlsls)
  studydata2[time - prevdur,] <- c(time - prevdur,casesls2,controlsls2)
  ###################################################################################################
  ### Stop when no more transmission 
  if (I1==0 & I2==0){
    break()
  }
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
  
  time <- time + 1
}

## saving the previous state
I1nvold2 <- I1nvold
I1vold2 <- I1vold
I2nvold2 <- I2nvold
I2vold2 <- I2vold
S1nvold2 <- S1nvold
S2nvold2 <- S2nvold
S0nvold2 <- S0nvold
S12nvold2 <- S12nvold
S1vold2 <- S1vold
S2vold2 <- S2vold
S12vold2 <- S12vold
## saving current state
I1nvold <- I1nv
I1vold <- I1v
I2nvold <- I2nv
I2vold <- I2v
S1nvold <- S1nv
S2nvold <- S2nv
S0nvold <- S0nv
S12nvold <- S12nv
S1vold <- S1v
S2vold <- S2v
S12vold <- S12v
## restoring state
I1nv <- I1nvold
I1v <- I1vold
I2nv <- I2nvold
I2v <- I2vold
S1nv <- S1nvold
S2nv <- S2nvold
S0nv <- S0nvold
S12nv <- S12nvold
S1v <- S1vold
S2v <- S2vold
S12v <- S12vold
## restoring previous state
I1nv <- I1nvold2
I1v <- I1vold2
I2nv <- I2nvold2
I2v <- I2vold2
S1nv <- S1nvold2
S2nv <- S2nvold2
S0nv <- S0nvold2
S12nv <- S12nvold2
S1v <- S1vold2
S2v <- S2vold2
S12v <- S12vold2

I1 <- sum(c(I1nv,I1v))
I2 <- sum(c(I2nv,I2v))

for (k in 1:1){
  ## virus 1
  p1inf <- pSE(r10,I1)
  inf11nv <- rbinom(1, S1nv, p1inf)
  inf21nv <- rbinom(1, S2nv, p1inf)
  inf01nv <- rbinom(1, S0nv, p1inf)
  inf121nv <- rbinom(1, S12nv, p1inf)
  
  inf11v <- sapply(S1v, function(sv) rbinom(1, sv, p1inf))
  inf11v2 <- sapply(seq_along(inf11v), function(k) ifelse(k <= (time - prevdur), 0, inf11v[k]))
  inf121v <- sapply(S12v, function(sv) rbinom(1, sv, p1inf))
  inf121v2 <- sapply(seq_along(inf121v), function(k) ifelse(k <= (time - prevdur), 0, inf11v[k]))
  
  rem1nv <- rbinom(1,I1nv,prem)
  rem1v <- sapply(I1v, function(iv) rbinom(1,iv, prem))
  ## Updating # susceptibles, because used again for virus 2
  ## unvaccinated
  S1nv <- S1nv - inf11nv
  S2nv <- S2nv - inf21nv
  S0nv <- S0nv - inf01nv
  S12nv <- S12nv - inf121nv
  ## vaccinated
  S1v <- S1v - inf11v
  S12v <- S12v - inf121v
  
  ## for plotting purposes ...
  I1ls <- c(I1ls,sum(inf11v) + sum(inf121v) + inf11nv + inf21nv + inf01nv + inf121nv)  
  ## virus 2
  p2inf <- pSE(r20,I2)
  inf12nv <- rbinom(1, S1nv, p2inf)
  inf22nv <- rbinom(1, S2nv, p2inf)
  inf02nv <- rbinom(1, S0nv, p2inf)
  inf122nv <- rbinom(1, S12nv, p2inf)
  
  inf22v <- sapply(S2v, function(sv) rbinom(1, sv, p2inf))
  inf22v2 <- sapply(seq_along(inf22v), function(k) ifelse(k <= (time - prevdur), 0, inf22v[k]))
  inf122v <- sapply(S12v, function(sv) rbinom(1, sv, p1inf))
  inf122v2 <- sapply(seq_along(inf122v), function(k) ifelse(k <= (time - prevdur), 0, inf11v[k]))
  
  rem2nv <- rbinom(1,I2nv,prem)
  rem2v <- sapply(I2v, function(iv) rbinom(1,iv, prem))
  ## Updating # susceptibles
  ## unvaccinated
  S1nv <- S1nv - inf12nv
  S2nv <- S2nv - inf22nv
  S0nv <- S0nv - inf02nv
  S12nv <- S12nv - inf122nv
  ## vaccinated
  S2v <- S2v - inf22v
  S12v <- S12v - inf122v
  
  ## for plotting purposes ...
  I2ls <- c(I2ls,sum(inf22v) + sum(inf122v) + inf12nv + inf22nv + inf02nv + inf122nv)  
  #########################################################################################
  ## updating other states: Infections with virus 1
  I1nv <- I1nv + inf11nv + inf21nv + inf01nv + inf121nv
  I1v <- I1v + inf11v + inf121v
  ## updating other states: Infections with virus 2
  I2nv <- I2nv + inf12nv + inf22nv + inf02nv + inf122nv
  I2v <- I2v + inf22v + inf122v
  ## updating states: Infections with either virus
  Rnv <- Rnv + rem1nv + rem2nv
  Rv <- Rv + rem1v + rem2v
  
  I1 <- sum(c(I1nv,I1v))
  I2 <- sum(c(I2nv,I2v))
  ###################################################################################################
  ### "Study" is conducted ##########################################################################
  ###################################################################################################
  infnv <- inf11nv + inf21nv + inf01nv + inf121nv + inf12nv + inf22nv + inf02nv + inf122nv
  infv <- inf11v + inf22v + inf121v + inf122v
  infv2 <- inf11v2 + inf22v2 + inf121v2 + inf122v2
  
  casesls <- c(infnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(infv[vaccdelim[x] : vaccdelim[x + 1]])))
  casesls2 <- c(infnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(infv2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  controlsvsum <- S1v + S2v + S0v + S12v + Rv
  controlsvsum2 <- sapply(seq_along(controlsvsum), function(k) ifelse(k <= (time - prevdur), 0, controlsvsum[k]))
  
  controlsnv <- S1nv + S2nv + S0nv + S12nv + Rnv
  
  controlsls <- c(controlsnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum[vaccdelim[x] : vaccdelim[x + 1]])))
  controlsls2 <- c(controlsnv,sapply(seq_along(vaccdelim[-1]), function(x) sum(controlsvsum2[vaccdelim[x] : vaccdelim[x + 1]])))
  
  studydata[time - prevdur,] <- c(time - prevdur,casesls,controlsls)
  studydata2[time - prevdur,] <- c(time - prevdur,casesls2,controlsls2)
  ###################################################################################################
  ### Stop when no more transmission 
  if (I1==0 & I2==0){
    break()
  }
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
  
  time <- time + 1
}

