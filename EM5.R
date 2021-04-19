
Estep <- function(d, dsub, oldp, indHs, BmPlas, BmTrts) {
  ns <- dim(dsub)[1]
  nhs <- length(oldp$h)
  L <- ncol(BmPlas)

  G <- matrix(exp(oldp$alpvec[1]*dsub$x), nrow=ns, ncol=nhs, byrow=FALSE)
  A0n <- matrix(dsub$A0, nrow=ns, ncol=nhs, byrow=FALSE)
  for(j in 1:L) G <- with(dsub, G*exp(oldp$alpvec[j+1]*((1-A0n)*matrix(BmPlas[,j],nrow=ns, ncol=nhs, byrow=TRUE)+A0n*matrix(BmTrts[,j],nrow=ns, ncol=nhs, byrow=TRUE))))     
    
  dsub$a <- with(dsub, 1/oldp$theta+Ni)
  dsub$b <- 1/(1/oldp$theta+ (indHs*G)%*%oldp$h)

  dsub$Exi <- with(dsub, a*b)
  dsub$Elogxi <- with(dsub, digamma(a)+log(b)) 
  dsub$Exi2 <- with(dsub, a*b^2+Exi^2) 
  dsub$Exilogxi <- with(dsub, a*b*digamma(a+1)+a*b*log(b)) 
  dsub$Elogxi2 <- with(dsub, trigamma(a+1)+Elogxi^2) 
  # cat("Estep 1"); gc()
  rm(list=ls()[ls() %nin% c("dsub")]) 
  #cat("Estep 2");  gc()

  return(dsub)
}                                              
  
######################################################################################################

Mstep <- function(d, dsub, oldp, indHs, BmPlas, BmTrts, indHNs, indSNs, tranf='iden') {
     Ns <- dim(d)[1]
     ns <- dim(dsub)[1]
     na <- length(oldp$alpvec)
     nhs <- length(oldp$h)
     dl1alp <- rep(NA, na)
     dl2alp <- matrix(NA, nrow=na, ncol=na)
 
     L <- ncol(BmPlas)
     G <- matrix(1, ns, nhs) 
     tmpN <- tmpn <- vector("list", na)
 
     if(tranf=='iden') {
        temptheta <- 1/oldp$theta
        newtheta <- temptheta - ((1+log(temptheta)-digamma(temptheta))+ with(dsub, mean(Elogxi-Exi)))/(1/temptheta-trigamma(temptheta))
        newtheta <- 1/newtheta
     }
 
     if(tranf=='log') {
        tt <- -log(oldp$theta)
        newtheta <- tt - ((exp(tt)*(tt+1)-digamma(exp(tt))*exp(tt))+ exp(tt)*with(dsub, mean(Elogxi-Exi)))/(exp(tt)*(tt+2)-trigamma(exp(tt))*exp(2*tt)-digamma(exp(tt))*exp(tt)+exp(tt)*with(dsub, mean(Elogxi-Exi)))
        newtheta <- exp(-newtheta)
     }

     A0N <- matrix(d$A0, nrow=Ns, ncol=nhs, byrow=FALSE)
     A0n <- matrix(dsub$A0, nrow=ns, ncol=nhs, byrow=FALSE)
     tmpN[[1]] <- matrix(d$x, nrow=Ns, ncol=nhs, byrow=FALSE)
     tmpn[[1]] <- matrix(dsub$x,ns,nhs,byrow=FALSE)
     for(j in 2:na) {
        tmpN[[j]] <- (1-A0N)*matrix(BmPlas[,j-1],Ns,nhs,byrow=TRUE)+A0N*matrix(BmTrts[,j-1],Ns,nhs,byrow=TRUE)
        tmpn[[j]] <- (1-A0n)*matrix(BmPlas[,j-1],ns,nhs,byrow=TRUE)+A0n*matrix(BmTrts[,j-1],ns,nhs,byrow=TRUE)
     }
    for(j in 1:na) G <- G*exp(oldp$alpvec[j]*tmpn[[j]])
    newh <- colSums(indHNs*indSNs)/colSums(indHs*G*matrix(dsub$Exi, nrow=ns, ncol=nhs, byrow=FALSE))
  
    for(j1 in 1:na) {
        dl1alp[j1] <- sum(indSNs*indHNs*tmpN[[j1]])-sum(colSums(indSNs*indHNs)*colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G*tmpn[[j1]])/colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G))
        for(j2 in j1:na) {
          dl2alp[j1,j2] <- dl2alp[j2,j1] <- -sum(colSums(indSNs*indHNs)*colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G*tmpn[[j1]]*tmpn[[j2]])/colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G))+
                            sum(colSums(indSNs*indHNs)*colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G*tmpn[[j1]])*colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G*tmpn[[j2]])/
                            (colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G))^2)   
      }
     }
     newalpvec <- oldp$alpvec - solve(dl2alp)%*%dl1alp
     newpara <- list(h=newh, alpvec=newalpvec, theta=newtheta)
     rm(list=ls()[ls() %nin% c("newpara")]) 
     return(newpara)
}

### This function is a modofied version of EM5.R file to incorporate staggered switching time

Estep.swit <- function(d, dsub, oldp, indHs, BmMats) {
 # revision of Estep() allowing different switching times across the patients 

  ns <- dim(dsub)[1]
  nhs <- length(oldp$h)
  L <- ncol(BmMat[[1]])
  mat <- matrix(NA, ns, nhs)
  
  G <- matrix(exp(oldp$alpvec[1]*dsub$x), nrow=ns, ncol=nhs, byrow=FALSE)
  for(j in 2:(L+1)) {
     for(i in 1:ns) mat[i,] <- BmMats[[i]][,j-1]     
     G <- with(dsub, G*exp(oldp$alpvec[j]*mat))     
  }
    
  dsub$a <- with(dsub, 1/oldp$theta+Ni)
  dsub$b <- 1/(1/oldp$theta+ (indHs*G)%*%oldp$h)

  dsub$Exi <- with(dsub, a*b)
  dsub$Elogxi <- with(dsub, digamma(a)+log(b)) 
  dsub$Exi2 <- with(dsub, a*b^2+Exi^2) 
  dsub$Exilogxi <- with(dsub, a*b*digamma(a+1)+a*b*log(b)) 
  dsub$Elogxi2 <- with(dsub, trigamma(a+1)+Elogxi^2) 
  # cat("Estep 1"); gc()
  rm(list=ls()[ls() %nin% c("dsub")]) 
  #cat("Estep 2");  gc()

  return(dsub)
}                                              
  
######################################################################################################

Mstep.swit <- function(d, dsub, oldp, indHs, BmMats, indHNs, indSNs, tranf='iden') {
 # revision of Mstep() allowing different switching times across the patients 

     Ns <- dim(d)[1]
     ns <- dim(dsub)[1]
     na <- length(oldp$alpvec)
     nhs <- length(oldp$h)
     dl1alp <- rep(NA, na)
     dl2alp <- matrix(NA, nrow=na, ncol=na)
 
     L <- ncol(BmMat[[1]])
     G <- matrix(1, ns, nhs) 
     matN <- matrix(NA, Ns, nhs)
     mat <- matrix(NA, ns, nhs)
     tmpN <- tmpn <- vector("list", na)
 
     if(tranf=='iden') {
        temptheta <- 1/oldp$theta
        newtheta <- temptheta - ((1+log(temptheta)-digamma(temptheta))+ with(dsub, mean(Elogxi-Exi)))/(1/temptheta-trigamma(temptheta))
        newtheta <- 1/newtheta
     }
 
     if(tranf=='log') {
        tt <- -log(oldp$theta)
        newtheta <- tt - ((exp(tt)*(tt+1)-digamma(exp(tt))*exp(tt))+ exp(tt)*with(dsub, mean(Elogxi-Exi)))/(exp(tt)*(tt+2)-trigamma(exp(tt))*exp(2*tt)-digamma(exp(tt))*exp(tt)+exp(tt)*with(dsub, mean(Elogxi-Exi)))
        newtheta <- exp(-newtheta)
     }

     tmpN[[1]] <- matrix(d$x, nrow=Ns, ncol=nhs, byrow=FALSE)
     tmpn[[1]] <- matrix(dsub$x,ns,nhs,byrow=FALSE)
     for(j in 2:na) {
        for(i in 1:Ns) matN[i,] <- BmMats[[d$ID[i]]][,j-1] #need better program, require d$id has value 1:n as dsub patient order
        for(i in 1:ns) mat[i,] <- BmMats[[i]][,j-1] 
        tmpN[[j]] <- matN
        tmpn[[j]] <- mat
     }
    for(j in 1:na) G <- G*exp(oldp$alpvec[j]*tmpn[[j]])
    newh <- colSums(indHNs*indSNs)/colSums(indHs*G*matrix(dsub$Exi, nrow=ns, ncol=nhs, byrow=FALSE))
  
    for(j1 in 1:na) {
        dl1alp[j1] <- sum(indSNs*indHNs*tmpN[[j1]])-sum(colSums(indSNs*indHNs)*colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G*tmpn[[j1]])/colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G))
        for(j2 in j1:na) {
          dl2alp[j1,j2] <- dl2alp[j2,j1] <- -sum(colSums(indSNs*indHNs)*colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G*tmpn[[j1]]*tmpn[[j2]])/colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G))+
                            sum(colSums(indSNs*indHNs)*colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G*tmpn[[j1]])*colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G*tmpn[[j2]])/
                            (colSums(indHs*matrix(dsub$Exi,ns,nhs,byrow=FALSE)*G))^2)   
      }
     }
     newalpvec <- oldp$alpvec - solve(dl2alp)%*%dl1alp
     newpara <- list(h=newh, alpvec=newalpvec, theta=newtheta)
     rm(list=ls()[ls() %nin% c("newpara")]) 
     return(newpara)
}

