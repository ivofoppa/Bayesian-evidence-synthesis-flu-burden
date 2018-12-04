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

a <- 25

R <- function(a){
  1 - (1-R0)*exp(-rho*a)
}

S <- function(a){
  R0*exp(-a*lambda) + rho*(1 - R0)/(lambda - rho) * (exp(-rho*a) - exp(-lambda*a))
}

not_inf <- function(a){
  1 - R(a) + R0*exp(-a*lambda) + rho*(1 - R0)/(lambda - rho) * (exp(-rho*a) - exp(-lambda*a))
}

not_inf_prime <- function(a){
   -rho*(1-R0)*exp(-rho*a) - lambda*R0*exp(-a*lambda) - rho*(1 - R0)/(lambda - rho) * (-rho*exp(-rho*a) + lambda*exp(-lambda*a))
}

- not_inf_prime(a)/not_inf(a)

### VE, based on "estimated inst. attack rate" at age a
1 - S(a)/not_inf(a)/(exp(-lambda*a))

VEls <- sapply(1:200,function(a) 1 - S(a)/not_inf(a)/(exp(-lambda*a)))

lambda * (R0*exp(-rho*a) + rho*(1 - R0)/(lambda - rho) * (exp(-rho*a) - exp(-lambda*a)))/
  (1 - R(a) + rho*(1 - R0)/(lambda - rho) * (exp(-rho*a) - exp(-lambda*a)))











