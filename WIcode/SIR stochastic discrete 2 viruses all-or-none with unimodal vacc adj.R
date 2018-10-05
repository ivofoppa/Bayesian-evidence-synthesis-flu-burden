## Discrete stochastic transmission model with unimdal vaccination uptake, two viruses, all-or-none protection

r10 <- 1.6
r20 <- 1.8
delta <- 1/4 ## infectious period

seas <- 500 ## Number of days in epidemic
### Generating the vaccination rate over time
vacc1 <- 0.47 ## cumulative vaccination coverage

vs1 <- .1 ## proportion of subject who, if vacc. remain susc. to virus 1
vs2 <- .4 ## proportion of subject who, if vacc. remain susc. to virus 2
vs0 <- .3 ## proportion of subject who, if vacc. remain susc. to neither virus
vs12 <- 1 - vs1 - vs2 - vs0 ## proportion of subject who, if vacc. remain susc. to both viruses

# ## scaling the second virus to about virus 1 transmissability
# fac1 <- vacc1 *vs1 + (1 - vacc1)
# fac2 <- vacc1 *vs2 + (1 - vacc1)
# r20 <- r10 * fac2/fac1

prevdur <- 50 ## number of days of vaccination campaigh before transmission
vdur <- prevdur + 200 ### duration of vaccination
vaccint <- 60
ccratio <- Inf ## control-case ratio
###################################################################################################
##  Seed infections ###############################################################################
###################################################################################################
inf1num <- 90 ## Number of seed infections
inf2num <- 10 ## Number of seed infections
###################################################################################################
## Setting initial values etc.
source('initialize 2-virus model.R') ## run initialization code

## Running model
source('model run.R') ## run initialization code

## ... and setting-up data
source('data set prep.R') ## run initialization code

cond_logist <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset,method = 'approximate')
summary(cond_logist)

cond_logist2 <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset2,method = 'approximate')
summary(cond_logist2)

logist <- glm(cbind(cases,controls) ~ vacc, data = dataset3,family = binomial(link = 'logit'))
summary(logist)
######################################################################################################
######################################################################################################
plot(VEls)

plot(I1ls, ylim = c(0,max(c(I1ls,I2ls))),type = 'l',col = 'blue')
lines(I2ls,col = 'red')

filepath <- paste0('C:/Users/VOR1/Documents/GitHub/Waning-Immunity-artefact/WIwriteup/WIplots/simul_adj_',prevdur,'_',vdur,'.RData')
save(dataset,dataset2,studydata,studydata2,file = filepath)

