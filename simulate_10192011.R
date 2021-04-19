
Solve.surv.C <- function(told, p, tol=1.0e-5, maxiter=200) { 

  absdiff <- 1
  iter <- 0
  if(p$A0==1) {
    cont <- -(1-p$Rmax)*p$T0-p$Rmax*p$T50*log(p$T50+p$T0)+log(p$U)*exp(-p$gamma*p$xs)/p$lambda0/p$xis
    while(absdiff>tol & iter<maxiter) {
      d1 <- (1-p$Rmax)*told+p$Rmax*p$T50*log(p$T50+told)+cont
      d2 <- 1-p$Rmax+p$Rmax*p$T50/(p$T50+told)
      tnew <- told-d1/d2
      absdiff <- abs(tnew-told)
      told <- tnew
      iter <- iter+1
    } 
  }

  if(p$A0==0) {
    if(p$T0<p$c0) cont <- -(1-p$Rmax)*p$c0-p$Rmax*p$T50*log(p$T50)+log(p$U)*exp(-p$gamma*p$xs)/p$lambda0/p$xis+p$c0-p$T0
    if(p$T0>=p$c0) cont <- -(1-p$Rmax)*p$T0-p$Rmax*p$T50*log(p$T50+p$T0-p$c0)+log(p$U)*exp(-p$gamma*p$xs)/p$lambda0/p$xis
    while(absdiff>tol & iter<maxiter) {
      d1 <- (1-p$Rmax)*told+p$Rmax*p$T50*log(p$T50+told-p$c0)+cont
      d2 <- 1-p$Rmax+p$Rmax*p$T50/(p$T50+told-p$c0)
      tnew <- told-d1/d2
      absdiff <- abs(tnew-told)
      told <- tnew
      iter <- iter+1
      } 
  }
  return(tnew)
}

############################################################################

Solve.surv.D <- function(told, p, tol=1.0e-5, maxiter=200) { 

  absdiff <- 1
  iter <- 0
  if(p$A0==1) {
    cont <- -p$T0+0.5*p$alpha0/p$alpha1*exp(-p$alpha1*p$T0^2)+log(p$U)*exp(-p$gamma*p$xs)/p$lambda0/p$xis
    while(absdiff>tol & iter<maxiter) {
      d1 <- told-0.5*p$alpha0/p$alpha1*exp(-p$alpha1*told^2)+cont
      d2 <- 1+p$alpha0*told*exp(-p$alpha1*told^2)
      tnew <- told-d1/d2
      absdiff <- abs(tnew-told)
      told <- tnew
      iter <- iter+1
    } 
  }

  if(p$A0==0) {
    if(p$T0<p$c0) cont <- -p$T0+0.5*p$alpha0/p$alpha1+log(p$U)*exp(-p$gamma*p$xs)/p$lambda0/p$xis
    if(p$T0>=p$c0) cont <- -p$T0+0.5*p$alpha0/p$alpha1*exp(-p$alpha1*(p$T0-p$c0)^2)+log(p$U)*exp(-p$gamma*p$xs)/p$lambda0/p$xis
    while(absdiff>tol & iter<maxiter) {
      d1 <- told-0.5*p$alpha0/p$alpha1*exp(-p$alpha1*(told-p$c0)^2)+cont
      d2 <- 1+p$alpha0*(told-p$c0)*exp(-p$alpha1*(told-p$c0)^2)
      tnew <- told-d1/d2
      absdiff <- abs(tnew-told)
      told <- tnew
      iter <- iter+1
      } 
  }
  return(tnew)
}

#######################################################################################

simulate.Const <- function(npla, ntrt, c0, t0, theta, phiv, gamma, lambda0, CentL, tau)
# simulation for scenario (a) and (b) with constant treatment effect
{

  n <- npla+ntrt
  A0 <-  c(rep(0, npla), rep(1,ntrt))
  c0t <- c0*(1-A0)
  x <- rpois(n, 2)
  x <- ifelse(x<=5, x, 5)
  xi <- rgamma(n, shape=1/theta, rate=1/theta)
  Xc <- cbind(rep(1,n), A0, x)
  Cent <- as.vector(-log(runif(n))*exp(-Xc%*%phiv)) # censoring time
  Cent <- ifelse(Cent<tau, ifelse(Cent<CentL, CentL, Cent), tau) # censoring time truncated by (CentL, tau) 

  ID <- Y <- Event <- X <- Seq <- At <- A0t <- Centt <- Ni <- NULL

  for (i in 1:n) { 
    tempY=NULL; lastY=0; lastA=A0[i]; newA=A0[i];
    while(lastY<Cent[i]) {
      # simulate the survival time
      newY <- (-log(runif(1))/(lambda0*xi[i])*exp(-x[i]*gamma-t0*lastA))+lastY
      if (newY>c0t[i] & lastA==0) {
          newA=1
          newY <- (newY-lastY-(c0t[i]-lastY))*exp(-t0)+c0t[i]
      }
      if (newY<Cent[i]) tempY <- c(tempY, newY) 
      lastA=newA
      lastY <- newY
      At <- c(At, lastA)
    }
    tempY <- c(tempY, Cent[i])
    ni <- length(tempY) 
    ID <- c(ID, rep(i, ni))
    Y <- c(Y, tempY)
    Seq <- c(Seq, 1:ni)
    Event <- c(Event, rep(1, ni-1), 0)
    X <- c(X, rep(x[i], ni))
    A0t <- c(A0t, rep(A0[i],ni))
    Ni <- c(Ni, rep(ni-1, ni))
    Centt <- c(Centt, rep(Cent[i], ni))
  }

  dat <- data.frame(ID=ID, Survt=Y, Trt=At, A0=A0t, x=X, Seq=Seq, Cent=Centt, CenI=1-Event, event=Event, Ni=Ni)
  d1 <- subset(dat, Seq==1, select=c("ID", "A0", "x", "Cent", "Ni")) # one row per subject
  rm(list=ls()[ls() %nin% c("dat", "d1")])

  return(list(duniq=d1, dwhole=dat))
}

###############################################################################################
simulate.C <- function(npla, ntrt, c0, T50, Rmax, theta, phiv, gamma, lambda0, CentL, tau)
# simulation for scenario (c) with treatment function of log(1-Rmax*t/(T50+t))
{

  n <- npla+ntrt
  A0 <-  c(rep(0, npla), rep(1,ntrt))
  c0t <- c0*(1-A0)
  x <- rpois(n, 2)
  x <- ifelse(x<=5, x, 5)
  xi <- rgamma(n, shape=1/theta, rate=1/theta)
  Xc <- cbind(rep(1,n), A0, x)
  Cent <- as.vector(-log(runif(n))*exp(-Xc%*%phiv)) # censoring time
  Cent <- ifelse(Cent<tau, ifelse(Cent<CentL, CentL, Cent), tau) # censoring time truncated by (CentL, tau) 

  ID <- Y <- Event <- X <- Seq <- At <- A0t <- Centt <- Ni <- NULL

  for (i in 1:n) { 
    tempY=NULL; lastY=0; lastA=A0[i]; newA=A0[i];
    while(lastY<Cent[i]) {
      # simulate the survival time
      U <- runif(1)
      if (newA==0) newY <- -log(U)*exp(-gamma*x[i])/xi[i]/lambda0+lastY
      if (newA==1) newY <- Solve.surv.C(told=lastY, p=list(Rmax=Rmax, T50=T50, T0=lastY, A0=A0[i], c0=c0t[i], xs=x[i], gamma=gamma, xis=xi[i], lambda0=lambda0, U=U), tol=1.0e-5, maxiter=200) 
       # BBsolve(par=t, fn=equfun, Rmax=Rmax, T50=T50, T0=lastY, xs=x[k], gamma=gamma, xis=xi[k], lambda0=lambda0, C=lastY, U=U)
      if (newY>c0t[i] & lastA==0) {
          newA=1
          newY <- Solve.surv.C(told=newY, p=list(Rmax=Rmax, T50=T50, T0=lastY, A0=A0[i], c0=c0t[i], xs=x[i], gamma=gamma, xis=xi[i], lambda0=lambda0, U=U), tol=1.0e-5, maxiter=200)
      }
      # cat("sub=", i, "newY=", newY, "\n");
      if (newY<Cent[i]) tempY <- c(tempY, newY) 
      lastA=newA
      lastY <- newY
      At <- c(At, lastA)
    }
    tempY <- c(tempY, Cent[i])
    ni <- length(tempY) 
    ID <- c(ID, rep(i, ni))
    Y <- c(Y, tempY)
    Seq <- c(Seq, 1:ni)
    Event <- c(Event, rep(1, ni-1), 0)
    X <- c(X, rep(x[i], ni))
    A0t <- c(A0t, rep(A0[i],ni))
    Ni <- c(Ni, rep(ni-1, ni))
    Centt <- c(Centt, rep(Cent[i], ni))
  }

  dat <- data.frame(ID=ID, Survt=Y, Trt=At, A0=A0t, x=X, Seq=Seq, Cent=Centt, CenI=1-Event, event=Event, Ni=Ni)
  d1 <- subset(dat, Seq==1, select=c("ID", "A0", "x", "Cent", "Ni")) # one row per subject
  rm(list=ls()[ls() %nin% c("dat", "d1")])
  return(list(duniq=d1, dwhole=dat))
}

#############################################################################################

simulate.D <- function(npla, ntrt, c0, alpha0, alpha1, theta, phiv, gamma, lambda0, CentL, tau)
# simulation for scenario (c) with treatment function of log(1-Rmax*t/(T50+t))
{

  n <- npla+ntrt
  A0 <-  c(rep(0, npla), rep(1,ntrt))
  c0t <- c0*(1-A0)
  x <- rpois(n, 2)
  x <- ifelse(x<=5, x, 5)
  xi <- rgamma(n, shape=1/theta, rate=1/theta)
  Xc <- cbind(rep(1,n), A0, x)
  Cent <- as.vector(-log(runif(n))*exp(-Xc%*%phiv)) # censoring time
  Cent <- ifelse(Cent<tau, ifelse(Cent<CentL, CentL, Cent), tau) # censoring time truncated by (CentL, tau) 

  ID <- Y <- Event <- X <- Seq <- At <- A0t <- Centt <- Ni <- NULL

  for (i in 1:n) { 
    tempY=NULL; lastY=0; lastA=A0[i]; newA=A0[i];
    while(lastY<Cent[i]) {
      # simulate the survival time
      U <- runif(1)
      if (newA==0) newY <- -log(U)*exp(-gamma*x[i])/xi[i]/lambda0+lastY
      if (newA==1) newY <- Solve.surv.D(told=lastY, p=list(alpha0=alpha0, alpha1=alpha1, T0=lastY, A0=A0[i], c0=c0t[i], xs=x[i], gamma=gamma, xis=xi[i], lambda0=lambda0, U=U), tol=1.0e-5, maxiter=200) 
       # BBsolve(par=t, fn=equfun, Rmax=Rmax, T50=T50, T0=lastY, xs=x[k], gamma=gamma, xis=xi[k], lambda0=lambda0, C=lastY, U=U)
      if (newY>c0t[i] & lastA==0) {
          newA=1
          newY <- Solve.surv.D(told=newY, p=list(alpha0=alpha0, alpha1=alpha1, T0=lastY, A0=A0[i], c0=c0t[i], xs=x[i], gamma=gamma, xis=xi[i], lambda0=lambda0, U=U), tol=1.0e-5, maxiter=200)
      }
      # cat("sub=", i, "newY=", newY, "\n");
      if (newY<Cent[i]) tempY <- c(tempY, newY) 
      lastA=newA
      lastY <- newY
      At <- c(At, lastA)
    }
    tempY <- c(tempY, Cent[i])
    ni <- length(tempY) 
    ID <- c(ID, rep(i, ni))
    Y <- c(Y, tempY)
    Seq <- c(Seq, 1:ni)
    Event <- c(Event, rep(1, ni-1), 0)
    X <- c(X, rep(x[i], ni))
    A0t <- c(A0t, rep(A0[i],ni))
    Ni <- c(Ni, rep(ni-1, ni))
    Centt <- c(Centt, rep(Cent[i], ni))
  }

  dat <- data.frame(ID=ID, Survt=Y, Trt=At, A0=A0t, x=X, Seq=Seq, Cent=Centt, CenI=1-Event, event=Event, Ni=Ni)
  d1 <- subset(dat, Seq==1, select=c("ID", "A0", "x", "Cent", "Ni")) # one row per subject
  rm(list=ls()[ls() %nin% c("dat", "d1")])
  return(list(duniq=d1, dwhole=dat))
}

#############################################################################################

simulate.E <- function(npla, ntrt, alpha0, alpha1, theta, phiv, gamma, lambda0, CentL, tau)
# simulation for scenario (E) with treatment function of log(1-Rmax*t/(T50+t))
{

  n <- npla+ntrt
  A0 <-  c(rep(0, npla), rep(1,ntrt))
  c0t <- runif(n, 0, 2)*(1-A0)
  x <- rpois(n, 2)
  x <- ifelse(x<=5, x, 5)
  xi <- rgamma(n, shape=1/theta, rate=1/theta)
  Xc <- cbind(rep(1,n), A0, x)
  Cent <- as.vector(-log(runif(n))*exp(-Xc%*%phiv)) # censoring time
  Cent <- ifelse(Cent<tau, ifelse(Cent<CentL, CentL, Cent), tau) # censoring time truncated by (CentL, tau) 

  ID <- Y <- Event <- X <- Seq <- At <- A0t <- Centt <- Ni <- c0 <- NULL

  for (i in 1:n) { 
    tempY=NULL; lastY=0; lastA=A0[i]; newA=A0[i];
    while(lastY<Cent[i]) {
      # simulate the survival time
      U <- runif(1)
      if (newA==0) newY <- -log(U)*exp(-gamma*x[i])/xi[i]/lambda0+lastY
      if (newA==1) newY <- Solve.surv.D(told=lastY, p=list(alpha0=alpha0, alpha1=alpha1, T0=lastY, A0=A0[i], c0=c0t[i], xs=x[i], gamma=gamma, xis=xi[i], lambda0=lambda0, U=U), tol=1.0e-5, maxiter=200) 
       # BBsolve(par=t, fn=equfun, Rmax=Rmax, T50=T50, T0=lastY, xs=x[k], gamma=gamma, xis=xi[k], lambda0=lambda0, C=lastY, U=U)
      if (newY>c0t[i] & lastA==0) {
          newA=1
          newY <- Solve.surv.D(told=newY, p=list(alpha0=alpha0, alpha1=alpha1, T0=lastY, A0=A0[i], c0=c0t[i], xs=x[i], gamma=gamma, xis=xi[i], lambda0=lambda0, U=U), tol=1.0e-5, maxiter=200)
      }
      # cat("sub=", i, "newY=", newY, "\n");
      if (newY<Cent[i]) tempY <- c(tempY, newY) 
      lastA=newA
      lastY <- newY
      At <- c(At, lastA)
    }
    tempY <- c(tempY, Cent[i])
    ni <- length(tempY) 
    ID <- c(ID, rep(i, ni))
    Y <- c(Y, tempY)
    Seq <- c(Seq, 1:ni)
    Event <- c(Event, rep(1, ni-1), 0)
    X <- c(X, rep(x[i], ni))
    A0t <- c(A0t, rep(A0[i],ni))
    Ni <- c(Ni, rep(ni-1, ni))
    Centt <- c(Centt, rep(Cent[i], ni))
    c0 <- c(c0, rep(c0t[i], ni))
  }

  dat <- data.frame(ID=ID, Survt=Y, Trt=At, A0=A0t, x=X, Seq=Seq, Cent=Centt, CenI=1-Event, event=Event, Ni=Ni, c=c0)
  d1 <- subset(dat, Seq==1, select=c("ID", "A0", "x", "Cent", "Ni", "c")) # one row per subject
  rm(list=ls()[ls() %nin% c("dat", "d1")])
  return(list(duniq=d1, dwhole=dat))
}
