
Covest <- function(dsub, oldp, indHs, BmPlas, BmTrts, indHNs, indSNs) {
  ns <- dim(dsub)[1]
  nhs <- length(oldp$h)
  na <- length(oldp$alpvec)
  tmpn <- vector("list", na)
  tmpv <- matrix(NA, ns, na)
  d2Q <- D <- matrix(0, (na+nhs+1),(na+nhs+1))
  G <- matrix(1, ns, nhs) 
 
  tmpn[[1]] <- matrix(dsub$x,ns,nhs,byrow=FALSE)
  A0n <- matrix(dsub$A0, nrow=ns, ncol=nhs, byrow=FALSE)
  for(j in 2:na) tmpn[[j]] <- (1-A0n)*matrix(BmPlas[,j-1],ns,nhs,byrow=TRUE)+A0n*matrix(BmTrts[,j-1],ns,nhs,byrow=TRUE)
  for(j in 1:na) G <- G*exp(oldp$alpvec[j]*tmpn[[j]])
    
  for (j1 in 1:na) {
	  for (j2 in j1:na) {
	  	d2Q[j2, j1] <- d2Q[j1, j2] <- with(dsub, -sum(Exi*((indHs*G*tmpn[[j1]]*tmpn[[j2]])%*%oldp$h)))    
    }
    d2Q[(na+1):(na+nhs),j1] <- d2Q[j1,(na+1):(na+nhs)] <- with(dsub, -as.vector(Exi)%*%(G*indHs*tmpn[[j1]]))
  }
  d2Q[(na+1):(na+nhs),(na+1):(na+nhs)] <- -diag(colSums(indHNs*indSNs)/(oldp$h)^2)  
  d2Q[(na+nhs+1), (na+nhs+1)] <- ns*(-trigamma(1/oldp$theta)+oldp$theta)
     
  for (j in 1:na) tmpv[,j] <- (indHs*G*tmpn[[j]])%*%oldp$h
  D[1:na,1:na] <- with(dsub, t(as.vector(a*b^2)*tmpv)%*%tmpv)
  D[1:na,(na+1):(na+nhs)] <- with(dsub, t(as.vector(a*b^2)*tmpv)%*%(G*indHs))
  D[1:na,(na+nhs+1)] <- with(dsub, colSums(as.vector(a*b^2-Exilogxi+Exi*Elogxi)*tmpv))
  D[(na+1):(na+nhs),(na+1):(na+nhs)] <- with(dsub, t(as.vector(a*b^2)*G*indHs)%*%(G*indHs))
  D[(na+1):(na+nhs),(na+nhs+1)] <- with(dsub, colSums(as.vector(a*b^2-Exilogxi+Exi*Elogxi)*(G*indHs)))
  D[(na+nhs+1), (na+nhs+1)] <- with(dsub, sum(a*b^2+Elogxi2-2*Exilogxi-Elogxi^2+2*Exi*Elogxi))
  D[(na+1):(na+nhs+1),1:na] <- t(D[1:na,(na+1):(na+nhs+1)])
  D[(na+nhs):(na+nhs+1),(na+1):(na+nhs)] <- t(D[(na+1):(na+nhs), (na+nhs):(na+nhs+1)])
 
  Cov <- solve(-d2Q-D)
  
  rm(list=ls()[ls() %nin% c("Cov")]) 
  return(Cov)
}         



Covest.swit <- function(dsub, oldp, indHs, BmMats, indHNs, indSNs) {
 # revision of Covest() allowing different switching times across the patients 
  ns <- dim(dsub)[1]
  nhs <- length(oldp$h)
  na <- length(oldp$alpvec)
  tmpn <- vector("list", na)
  tmpv <- matrix(NA, ns, na)
  d2Q <- D <- matrix(0, (na+nhs+1),(na+nhs+1))
  G <- matrix(1, ns, nhs) 
  mat <- matrix(NA, ns, nhs)
  
  tmpn[[1]] <- matrix(dsub$x,ns,nhs,byrow=FALSE)
  for(j in 2:na) {
    for(i in 1:n) mat[i,] <- BmMats[[i]][,j-1]
    tmpn[[j]] <- mat
  }
  for(j in 1:na) G <- G*exp(oldp$alpvec[j]*tmpn[[j]])
    
  for (j1 in 1:na) {
          for (j2 in j1:na) {
                d2Q[j2, j1] <- d2Q[j1, j2] <- with(dsub, -sum(Exi*((indHs*G*tmpn[[j1]]*tmpn[[j2]])%*%oldp$h)))    
    }
    d2Q[(na+1):(na+nhs),j1] <- d2Q[j1,(na+1):(na+nhs)] <- with(dsub, -as.vector(Exi)%*%(G*indHs*tmpn[[j1]]))
  }
  d2Q[(na+1):(na+nhs),(na+1):(na+nhs)] <- -diag(colSums(indHNs*indSNs)/(oldp$h)^2)  
  d2Q[(na+nhs+1), (na+nhs+1)] <- ns*(-trigamma(1/oldp$theta)+oldp$theta)
     
  for (j in 1:na) tmpv[,j] <- (indHs*G*tmpn[[j]])%*%oldp$h
  D[1:na,1:na] <- with(dsub, t(as.vector(a*b^2)*tmpv)%*%tmpv)
  D[1:na,(na+1):(na+nhs)] <- with(dsub, t(as.vector(a*b^2)*tmpv)%*%(G*indHs))
  D[1:na,(na+nhs+1)] <- with(dsub, colSums(as.vector(a*b^2-Exilogxi+Exi*Elogxi)*tmpv))
  D[(na+1):(na+nhs),(na+1):(na+nhs)] <- with(dsub, t(as.vector(a*b^2)*G*indHs)%*%(G*indHs))
  D[(na+1):(na+nhs),(na+nhs+1)] <- with(dsub, colSums(as.vector(a*b^2-Exilogxi+Exi*Elogxi)*(G*indHs)))
  D[(na+nhs+1), (na+nhs+1)] <- with(dsub, sum(a*b^2+Elogxi2-2*Exilogxi-Elogxi^2+2*Exi*Elogxi))
  D[(na+1):(na+nhs+1),1:na] <- t(D[1:na,(na+1):(na+nhs+1)])
  D[(na+nhs):(na+nhs+1),(na+1):(na+nhs)] <- t(D[(na+1):(na+nhs), (na+nhs):(na+nhs+1)])
  
  Cov <- solve(-d2Q-D)
  
  rm(list=ls()[ls() %nin% c("Cov")]) 
  return(Cov)
}                                              
                                     
