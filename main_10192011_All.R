
rm(list = ls(all = TRUE))
source("EM5.R")
source("simulate_10192011.R")
source("Covest.R")

input <-  read.table("input5.csv", skip = 1, header = FALSE, nrows=1)
index <- as.integer(input[1])
seed <- as.integer(input[2])
t0 <- as.double(input[3])
M <- as.integer(input[4])
nsim <- as.integer(input[5])

library(survival)
require(splines)
library(Hmisc)

nsim <- 5
n = 400; n0=200; n1=200;
Kn=6; 
# M=3; # the order of spline, if 3, the spline is quadratic
set.seed(seed)

# set parameters
gamma =0.2; theta = 0.8; lambda0 = 1.5; tau=3; alpha0=-0.5; alpha1=0.75;
phiv = c(1, 1, 1); phiv =c(-2, -0.2, 0.2);
Rmax= 0.9; T50 = 0.6; 
c0 = 1; # length of core study

# save results
Npla <- Ntrt <- rep(NA, nsim) # number of average events in each group
Surv.est <- matrix(NA, nrow=nsim, ncol=2)
Spli.est <- matrix(NA, nrow=nsim, ncol=1+Kn+M+1)
Spli.sd <- Spli.est <- matrix(NA, nrow=nsim, ncol=1+Kn+M+1)
covbetav <- matrix(NA, nrow=nsim, ncol=(Kn+M)*(Kn+M))


for (m in 1:nsim) {
  print(Sys.time())
 # simulate data and saved them into simdat
  if(index==1 | index==2) tmp <- simulate.Const(npla=n0, ntrt=n1, c0=c0, t0=t0, theta=theta, phiv=phiv, gamma=gamma, lambda0=lambda0, CentL=1.2, tau=tau)
  if(index==3) tmp <- simulate.C(npla=n0, ntrt=n1, c0=c0, T50=T50, Rmax=Rmax, theta=theta, phiv=phiv, gamma=gamma, lambda0=lambda0, CentL=1.2, tau=tau)
  if(index==4) tmp <- simulate.D(npla=n0, ntrt=n1, c0=c0, alpha0=alpha0, alpha1=alpha1, theta=theta, phiv=phiv, gamma=gamma, lambda0=lambda0, CentL=1.2, tau=tau)
  if(index==5) tmp <- simulate.E(npla=n0, ntrt=n1, alpha0=alpha0, alpha1=alpha1, theta=theta, phiv=phiv, gamma=gamma, lambda0=lambda0, CentL=1.2, tau=tau)
  d1 <- tmp$duniq
  dat <- tmp$dwhole
  rm(tmp)
  Surv.est[m, ] <- coxph(Surv(Survt, event) ~ x+ A0+ frailty(ID, method='em'),dat)$coefficients

  Npla[m] <- nrow(dat[dat$A0==0,])/n0-1
  Ntrt[m] <- nrow(dat[dat$A0==1,])/n1-1 
  ht <- unique(sort(with(subset(dat, event==1), Survt)))
  nh <- length(ht) # nh is the number of unique event times
  N <- dim(dat)[1] # N is the total number of time including event time and censoring time for all subjects, each subject may have multiple time points
  indH <- (matrix(d1$Cent, n, nh, byrow=FALSE) >= matrix(ht, n, nh, byrow=TRUE)) # I(Y_i >= lambda(t)) n*nh matrix, one row per subject 
  indHN <- (matrix(dat$Cent, N, nh, byrow=FALSE) >= matrix(ht, N, nh, byrow=TRUE)) # I(Y_{ij} >= lambda(t)) N*nh matrix, multiple rows per subject 
  indSN <- (matrix(dat$Survt, N, nh, byrow=FALSE) == matrix(ht, N, nh, byrow=TRUE)) # I(EventTime_{ij} == lambda(t)) N*nh matrix, multiple rows per subject
  indSN[dat$CenI==1,] <- 0 # in case the censoring time is same as event time

  knots <- c(rep(0,M-1), seq(0,tau,tau/(Kn+1)),rep(tau,M-1))
  if(index %in% c(1,2,3,4)){
    BmPla <- splineDesign(knots, ht-c0, M, outer.ok=TRUE) # generate B-spline basis functions for each event time for control group, nh*9 matrix
    BmTrt <- splineDesign(knots, ht, M, outer.ok=TRUE) # generate B-spline basis functions for each event time for treatment group, nh*9 matrix
  }
  if(index %in% c(5)){
    BmMat <-  vector("list", n)
    for (i in 1:n) BmMat[[i]] <- splineDesign(knots, ht-d1$c[i], M, outer.ok=TRUE) # generate B-spline basis functions for each event time for ith patient, nh*9 matrix
  }
  
  oldh <- rep(1/nh, nh)  
  oldalpvec <- c(Surv.est[m,1], rep(Surv.est[m,2], Kn+M)) # include gamma 
  # oldalpvec <- c(0.2, rep(0, Kn+M)) # include gamma 
  oldtheta <- 1
  oldpara <- list(h=oldh, alpvec=oldalpvec, theta=oldtheta)
  
  epsilon <- 0.001
  maxiter <- 200
  absdiff <- 1
  iter <- 0
  diffvec <- rep(1, length(oldpara))
  
  while(absdiff>epsilon & iter<maxiter) {
    # cat("alpha=", oldpara$alpvec, "\n", "theta=", oldpara$theta, "\n", sep="\t") 
    if(index %in% c(1,2,3,4)) d1 <- Estep(d=dat, dsub=d1, oldp=oldpara, indHs=indH, BmPlas=BmPla, BmTrts=BmTrt)
    if(index %in% c(5)) d1 <- Estep.swit(d=dat, dsub=d1, oldp=oldpara, indHs=indH, BmMats=BmMat)
    if(index %in% c(1,3,4)) newpara <- Mstep(d=dat, dsub=d1, oldp=oldpara, indHs=indH, BmPlas=BmPla, BmTrts=BmTrt, indHNs=indHN, indSNs=indSN, tranf='log')
    if(index %in% c(2)) newpara <- Mstep(d=dat, dsub=d1, oldp=oldpara, indHs=indH, BmPlas=BmPla, BmTrts=BmTrt, indHNs=indHN, indSNs=indSN, tranf='iden')
    if(index %in% c(5)) newpara <- Mstep.swit(d=dat, dsub=d1, oldp=oldpara, indHs=indH, BmMats=BmMat, indHNs=indHN, indSNs=indSN, tranf='log')
    for (k in 1:length(oldpara)) diffvec[k] <- max(abs(oldpara[[k]]-newpara[[k]]))
    absdiff <- max(diffvec)
    oldpara <- newpara
    iter <- iter+1
  }
  if(index %in% c(1,2,3,4)) Covmat <- Covest(dsub=d1, oldp=oldpara, indHs=indH, BmPlas=BmPla, BmTrts=BmTrt, indHNs=indHN, indSNs=indSN)
  if(index %in% c(5)) Covmat <- Covest.swit(dsub=d1, oldp=oldpara, indHs=indH, BmMats=BmMat, indHNs=indHN, indSNs=indSN)
  var <- diag(Covmat)
  npara <- length(var)

  covbeta <- Covmat[2:(Kn+M+1), 2:(Kn+M+1)]
  covbetav[m,] <- as.vector(covbeta)

  cat("index=", index, "sim=", m, "npara=", npara, "\n", sep="\t") 
  cat("# of negative variance", sum(var<0), "which are", which(var<0), "\n", sep="\t") 
  Spli.est[m,] <- c(newpara$alpvec, newpara$theta)
  Spli.sd[m,] <- sqrt(abs(var[c(1:(1+Kn+M), npara)])) 

  sim.data <- data.frame(indexV=rep(index, nsim), Surv.est, Spli.est, Spli.sd, Npla, Ntrt) 
  cov.data <- data.frame(indexV=rep(index, nsim), covbetav) 
  write.table(sim.data, file="SimData.txt", sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(cov.data, file="SimCov.txt", sep="\t", row.names=FALSE, col.names=FALSE)
}

# time <- seq(0,3,by=0.01)
# plot(x=time, y=splineDesign(knots, time, M)%*%newpara$alpvec[-1], xlim=c(0,3), ylim=c(-2,2), ylab="Fitted B-Spline Curve", main="Fitted B-Spline Curve -- Scenario (C)", type="l")
# lines(x=time, y=log(1-Rmax*time/(T50+time)), col='red', lty=1)


  sim.data <- data.frame(indexV=rep(index, nsim), Surv.est, Spli.est, Spli.sd, Npla, Ntrt) 
  cov.data <- data.frame(indexV=rep(index, nsim), covbetav) 
  write.table(sim.data, file="SimData.txt", sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(cov.data, file="SimCov.txt", sep="\t", row.names=FALSE, col.names=FALSE)
