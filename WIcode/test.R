vacc <- .5

vdur <- 90 ### duration of vaccination
prevdur <- 50 ## vaccination before transmission

vrate <- -log(1 - vacc)/vdur

v <- c(rep(vrate,vdur),rep(0,seas - vdur))
vaccdelim <- c(seq(1,60,12),seas)
nvaccat <- length(vaccdelim) ## number of time-since-vacc categories

time <- 1

while (time <= prevdur){
  ###################################################################################################
  ###################################################################################################
  ### Vaccination
  Sv <- shift(Sv,1)
  svacc <- rbinom(1,Snv,1 - exp(-v[time]))
  Sv[1] <- svacc; Snv <- Snv - svacc
  
  Ev <- shift(Ev,1)
  evacc <- rbinom(1,Env,1 - exp(-v[time]))
  Ev[1] <- evacc; Env <- Env - evacc
  
  Iv <- shift(Iv,1)
  ivacc <- rbinom(1,Inv,1 - exp(-v[time]))
  Iv[1] <- ivacc; Inv <- Inv - ivacc
  
  time <- time + 1
}
