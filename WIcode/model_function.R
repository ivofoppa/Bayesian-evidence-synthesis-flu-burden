SIRvaccmod <- function(r10=1.5,r20=1.5,delta=0.25,seas=300,vacc1=0.47,vs=c(.2,.6,.1,.1),prevdur,vdur,vaccint,
                       ccratio=Inf,infnum1=90,infnum2=10){
  if(sum(vs)!=1){
    stop('The proportions (argument vs) have to add up to 1!')
  }
  vs1 <- vs[1] ## proportion of subject who, if vacc. remain susc. to virus 1
  vs2 <- vs[2] ## proportion of subject who, if vacc. remain susc. to virus 2
  vs0 <- vs[3] ## proportion of subject who, if vacc. remain susc. to neither virus
  vs12 <- vs[4] ## proportion of subject who, if vacc. remain susc. to both viruses
  ###################################################################################################
  ## Setting initial values etc.
  source('initialize 2-virus model.R') ## run initialization code
  
  ## Running model
  source('model run.R') ## run initialization code
  
  ## ... and setting-up data
  source('data set prep.R') ## run initialization code
  
  cond_logist <- clogit(case ~ sincevacc + strata(time),weights = count, data = dataset,method = 'approximate')
  modsum <- coef(summary(cond_logist))
  
  output_list <- list()
  output_list[['parameters']] <- c(r10=r10,r20=r20,delta=delta,seas=seas,vacc1=vacc1,vs=vs,prevdur=prevdur,vdur=vdur,vaccint=vaccint,
  ccratio=ccratio,infnum1=infnum1,infnum2=infnum2)
  output_list[['exp_coef']] <- modsum[,2]
  output_list[['pvalues']] <- modsum[,5]
  
  return(output_list)
}
