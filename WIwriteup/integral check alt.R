###################################################################################################
###  This code is designed to investigate which of the multiple solutions to the       ############
###  lem" from Farrington et al. 1992                                                  ############
###  Ivo M. Foppa, October 2018                                                        ############
###################################################################################################
Npop <- 2e+6 ### Source population
vacc <- .5   ### Vaccination rate
AR <- .2     ### Seasonal attack rate
seas <- 150  ### Duration of season in days
lambda <- -log(AR)/seas ### Daily infection rate (per person per day)
rho <- 0.01  ### Rate of vaccination protection loss
R0 <- .3     ### Initial vaccine failure (proportion)
###################################################################################################
###################################################################################################
###################################################################################################
R <- sapply(seq_along(times)[-1], function(d) R[d] <<- 1 - (1 - R[d-1])*exp(-rho));
Rprime <- diff(R)

surv_a <- function(a,u){
  Rprime[u]*exp(-lambda*(a - u))
}
###################################################################################################
###################################################################################################
F <- function(a) {
  (1-R0)*rho/(lambda - rho) * (exp(-rho*a) - exp(-lambda*a))
}

Fprime <- function(a) {
  (1-R0)*rho/(lambda - rho) * (-rho*exp(-rho*a) + lambda*exp(-lambda*a))
}

R <- function(a){
  1 - (1-R0)*exp(-rho*a)
}

Rprime <- function(a){
  rho*(1-R0)*exp(-rho*a)
}

S <- function(a){
  1 - R(a) + R0*exp(-lambda*a) + F(a)
}

Sprime <- function(a){
  -Rprime(a) -lambda*exp(-lambda*a) + Fprime(a)
}

Salt <- function(a){
  R0*exp(-lambda*a) + F(a)
}

Saltprime <- function(a){
  -lambda*exp(-lambda*a) + Fprime(a)
}

a <- 20
-Sprime(a)/S(a) ### Instantanous risk in vaccinated at "age" 20, according to Farrington
-Saltprime(a)/Salt(a) ### Instantanous risk in vaccinated at "age" 20, according to Farrington

lambda*(R0*exp(-lambda*a) + F(a))
F(0)











