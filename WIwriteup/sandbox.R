ls <- seq(0,200)
arr <- NULL

for (n in ls){
  for (m in ls){
    arr <- rbind(arr,c(n,m),deparse.level = 0)
  }
}

lambda1 <- 50; lambda2 <- 60

# lsa <- sapply(seq_along(arr[,1]),function(k) arr[k,1]*dpois(arr[k,1],lambda1)*dpois(arr[k,2],lambda2))

lsa <- sapply(seq_along(arr[,1]),function(k) 1/arr[k,2]*dpois(arr[k,1],lambda1)*dpois(arr[k,2],lambda2))


lsb <- sapply(seq_along(arr[,1]),function(k) arr[k,1]/arr[k,2]*dpois(arr[k,1],lambda1)*dpois(arr[k,2],lambda2))

delta <- 0.0001
lsc <- sapply(seq_along(arr[,1]),function(k) (arr[k,1] + delta)/(arr[k,2] + delta)*dpois(arr[k,1],lambda1)*dpois(arr[k,2],lambda2))

delind <- which(arr[,2]==0)
sum(lsa[-delind])
sum(lsc[-delind])

1/lambda2
1/(lambda2/(1 - exp(-lambda2)))

dpois(100,lambda)

Npop <- 1e+6

vacc <- .4
csfrac <- 0.001 ### case sampling fraction
cnfrac <- csfrac*.1

nsim <- 1000000
N1 <- rbinom(nsim,Npop,vacc)
N0 <- Npop - N1

lambda0 <- -log(0.85) ### 15% seasonal attack rate

xls <- seq(-3.5,3.5,length.out = 150)

lambdat0 <- dnorm(xls)
lambdat0 <- lambdat0 - min(lambdat0)
lambdat0 <- lambdat0/sum(lambdat0)*lambda0

casetimels0 <- NULL
N0 <- N0[100]
for (k in seq_along(lambdat0)){
  pinf <- 1 - exp(-lambdat0[k])
  ct <- rbinom(1,N0,pinf)
  casetimels0 <- c(casetimels0,ct)
  N0 <- N0 - ct
}

phi <- .6 ### VE
lambda1 <- lambda0*(1-phi)

cases1 <- rpois(nsim,N1*lambda1*csfrac)
cases0 <- rpois(nsim,N0*lambda0*csfrac)

controls1 <- rpois(nsim,N1*cnfrac)
controls0 <- rpois(nsim,N0*cnfrac)

mean(1 - cases1/cases0/(controls1/controls0))

t1ls <- sapply(seq_along(cases1),function(k) (N1[k] - cases1[k]) + sum(rexp(cases1[k]),lambda1)
inft0ls <- rexp(cases0,lambda0)
inft1ls[1:10]
