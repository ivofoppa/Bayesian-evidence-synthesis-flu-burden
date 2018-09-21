VE<-.6
RR<-1-VE #flu incidence rate ratio for vacc group

Nb<-c(0,0)
Nb[1]<-100000
Nb[2]<-20000 #number at baseline in high and low risk groups for each vaccination treatment (early and late); 
#must be equal for early and late by no unmeas confounding

r<-c(0,0)
r[1]<-.001
r[2]<- r[1] 
#nonflu incidence per day in two groups (low and high flu incidence). May be unequal among the groups and this changes results

f<-c(0,0)
f[1]<-.001
f[2]<-.004
#flu incidence per day in two risk groups

d<-c(0,0)
d[1]<-0
d[2]<-60
#day of vaccination for early and late groups, relative to flu season start

dt<-75 
#day on which incidence is compared between early and late vaccinees

elig=matrix(c(0,0,0,0),nrow=2,ncol=2)
a=elig
incflu = elig
incnon=elig

for (i in (1:2)){    #risk group 1=low, 2=hi
  for (j in (1:2))   #vaccination time 1=early, 2=late
  {elig[i,j] <-Nb[i]*exp(-d[j]*(f[i]+r[i])) #still eligible for inclusion in study at time of vaccination because no prior flu test
  a[i,j]<-elig[i,j]*exp(-(dt-d[j])*(f[i]*RR+r[i])) #at risk on day dt because no flu test between vax and day t
  incflu[i,j]<-a[i,j]*f[i]*RR # incident flu cases on day dt
  incnon[i,j]<-a[i,j]*r[i] # incident nonflu cases on day dt 
  }}

incnontot<-c(0,0)
incflutot<-c(0,0)
for (j in 1:2) {incflutot[j]<-incflu[1,j]+incflu[2,j]
incnontot[j]<-incnon[1,j]+incnon[2,j]}
OR<-incflutot[1]*incnontot[2]/incflutot[2]/incnontot[1] #comparing early vs. late vaccination; early worse = waning = OR>1
OR
