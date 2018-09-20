
vacc <- .4

campdur <- 50
vrate <- -log(1 - vacc)/campdur

v <- c(.10,rep(vrate,campdur),rep(0,durEpidemic - 51))
