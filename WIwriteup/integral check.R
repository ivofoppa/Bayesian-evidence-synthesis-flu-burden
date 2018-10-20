###################################################################################################
###  This code is designed to investigate which of the multiple solutions to the       ############
###  lem" from Farrington et al. 1992                                                  ############
###  Ivo M. Foppa, October 2018                                                        ############
###################################################################################################
Npop <- 2e+6 ### Source population
vacc <- .5   ### Vaccination rate
AR <- .2     ### Seasonal attack rate
seas <- 150  ### Duration of season in days
delta <- .01 ### size of the time step (in days)

time <- seq(0,seas,delta)

lambda <- -log(AR)/seas ### Daily infection rate (per person per day)
rho <- 0.01  ### Rate of vaccination protection loss
R0 <- .3     ### Initial vaccine failure (proportion)
###################################################################################################
###################################################################################################
Nv0 <- Npop*vacc   ### Numbers vaccinated
Nu0 <- Npop - Nv0  ### Numbers unvaccinated
###################################################################################################
###################################################################################################
Sv0 <- Nv0*(1 - R0) ### Number vaccinated and susceptible
Su0 <- Nu0          ### Number unvaccinated and susceptible
###################################################################################################
###################################################################################################
Sv <- Sv0           ### Initialize vaccinated susceptibles
Su <- Su0           ### Initialize unvaccinated susceptibles
R <- R0             ### Imnitialize vaccine failures over time
###################################################################################################
###################################################################################################
###################################################################################################
R <- sapply(seq_along(time)[-1], function(d) R[d] <<- 1 - (1 - R[d-1])*exp(-rho*delta));
Rprime <- diff(R)/delta

surv_a <- function(a,u){
  if (!(u %in% time) | !(a %in% time)){
    stop('The time argument have to be contained in the "time" object!')
  }
  Rpr <- Rprime[which(time==u)-1]
  Rpr*exp(-lambda*(a - u))
}

a <- 30
Fa <- sapply(time[2:which(time==a)],function(u) surv_a(a,u))
c(sum(Fa*delta),rho*(1 - R0)/(lambda - rho) * (exp(-rho*a) - exp(-lambda*a)))

Fals <- NULL

for (a in time[-1]){
  Fa <- sapply(time[2:which(time==a)],function(u) surv_a(a,u))
  Fals <- rbind(Fals,c(sum(Fa*delta),rho*(1 - R0)/(lambda - rho) * (exp(-rho*a) - exp(-lambda*a))),deparse.level = 0)
  if (a %in% 1:seas){
    cat(paste0('Up to age ',a,' done!\n'))
  }
}

Faprime <- diff(Fals[,1])/delta

a <- 30

rho*(1 - R0)/(lambda - rho) * (-rho*exp(-rho*a) + lambda*exp(-lambda*a))
Faprime[which(time==a)-1]
###################################################################################################
###################################################################################################
Sprimea <- -Rprime[which(time==a) - 1] - lambda*R0*exp(-a*lambda) - rho*(1 - R0)/(lambda - rho) * (-rho*exp(-rho*a) + lambda*exp(-lambda*a))

Sa <- 1 - R[which(time==a)] + R0*exp(-a*lambda) + rho*(1 - R0)/(lambda - rho) * (exp(-rho*a) - exp(-lambda*a))

-Sprimea/Sa











