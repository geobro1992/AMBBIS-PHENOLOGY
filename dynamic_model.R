######################################################################
# This file contains the functions required to fit the dynamic models
######################################################################
# The main functions are
# dynamic_ll_func - the likelihood function
# start_val_func - choose random starting values
# fit_dynamic - wrapper for fitting the model from  multiple starts
# fit_func - fit the dynamic model and organise the output
######################################################################

# Link functions
expit <- function(x){ exp(x) / (1 + exp(x)) }
logit <- function(x){ log(x / (1 - x)) }


# Function for deriving array of arrival proportions in the stopover formulation
bettafunc <- function(j, betta.est, a1, phi.est) {
  temp <- a1[, , j]
  for (b in 1:(j-1)) {
    temp <- temp + betta.est[, , b] * phi.est ^ length(b:(j-1))
  }	
  return(temp)
}	

# Likelihood function		
dynamic_ll_func <- function(parm) {	
  
  par.index <- 0
  switch(M,
         "1" = { # univoltine species
           # productivity can be constant, time varying, or varying with one or more covariates (additively)
           switch(rho.m,
                  const = {
                    rho.est <- matrix(exp(parm[par.index + 1]),nrow = nS, ncol = nY - 1); par.index <- par.index + 1},
                  indyear = {
                    rho.est <- matrix(exp(parm[(par.index + 1):(par.index + nY - 1)]), nrow = nS, ncol = nY - 1, byrow = TRUE); par.index <- par.index + nY - 1},
                  l1cov = {
                    rho.est <- exp(parm[par.index + 1] + parm[par.index + 2] * rho1_cov1); par.index <- par.index + 2},
                  l2cov = {
                    rho.est <- exp(parm[par.index + 1] + parm[par.index + 2] * rho1_cov1 + parm[par.index + 3] * rho1_cov2); par.index <- par.index + 3})
           
           # mu can be constant, time varying, or varying with one or more covariates (additively)
           switch(mu.m,
                  const = {
                    mu.est <- matrix(exp(parm[par.index + 1]),nrow = nS, ncol = nY); par.index <- par.index + 1},
                  indyear = {
                    mu.est <- matrix(exp(parm[(par.index + 1):(par.index + nY)]), nrow = nS, ncol = nY, byrow = TRUE); par.index <- par.index + nY},
                  l1cov = {
                    mu.est <- exp(parm[par.index + 1] + parm[par.index + 2] * mu1_cov1); par.index <- par.index + 2},
                  l2cov = {
                    mu.est <- exp(parm[par.index + 1] + parm[par.index + 2] * mu1_cov1 + parm[par.index + 3] * mu1_cov2); par.index <- par.index + 3})
           
           # sigma can be constant, time varying, or varying with one or more covariates (additively)		
           switch(sigma.m,
                  const = {
                    sigma.est <- matrix(exp(parm[par.index+1]),nrow=nS,ncol=nY); par.index <- par.index+1},
                  indyear = {
                    sigma.est <-  matrix(exp(parm[(par.index+1):(par.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE); par.index <- par.index+nY},
                  l1cov = {
                    sigma.est <- exp(parm[par.index+1] +  parm[par.index+2]*sigma1_cov1); par.index <- par.index+2},
                  l2cov = {
                    sigma.est <- exp(parm[par.index+1] +  parm[par.index+2]*sigma1_cov1 + parm[par.index+3]*sigma1_cov2); par.index <- par.index+3})
           
           switch(atype,
                  "N"={
                    a <- array(NA,dim=c(nY,nT,nS))
                    for(kyear in 1:nY){
                      a[kyear,,] <- matrix(dnorm(rep(1:nT,nS),rep(mu.est[,kyear],each=nT),sigma.est[,kyear]),nrow=nT,ncol=nS)}},
                  "SO"={
                    # survival probability can be constant, time varying, or varying with one or more covariates (additively)
                    switch(phi.m,
                           const={
                             phi.est <-  matrix(expit(parm[par.index+1]),nrow=nS,ncol=nY); par.index <- par.index + 1},
                           indyear={
                             phi.est <- matrix(expit(parm[(par.index+1):(par.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE); par.index <- par.index + nY},
                           l1cov={
                             phi.est <- expit(parm[par.index + 1] + parm[par.index + 2]*phi1_cov1); par.index <- par.index + 2},
                           l2cov={
                             phi.est <- expit(parm[par.index + 1] + parm[par.index + 2]*phi1_cov1 + parm[par.index + 3]*phi1_cov2); par.index <- par.index + 3})
                    
                    # Array for betta, describing the proportions of N arriving
                    betta.est <- array(c(pnorm(1,mean=mu.est,sd=sigma.est),pnorm(rep(2:(nT-1),each=nS*nY),mean=rep(mu.est,length(2:(nT-1))),sd=rep(sigma.est,length(2:(nT-1))))-pnorm(rep(1:(nT-2),each=nS*nY),mean=rep(mu.est,length(1:(nT-2))),sd=rep(sigma.est,length(1:(nT-2)))),1-pnorm(nT-1,mean=mu.est,sd=sigma.est)),c(nS,nY,nT))
                    a1 <- betta.est
                    a <- array(cbind(c(betta.est[,,1]),sapply(2:nT,FUN=bettafunc,betta.est,a1,phi.est)),c(nS,nY,nT))
                    a <- aperm(a,c(2,3,1))})
           
           a[is.na(y)] <- NA
           # Estimate N_{i,1} (equation 6) for use in the concentrated likelihood
           ysum <- apply(y,3,sum,na.rm=TRUE)
           asum <- apply(a[1,,],2,sum,na.rm=TRUE) + apply(a[2,,],2,sum,na.rm=TRUE)*rho.est[,1]
           for(kyear in 3:nY){
             asum <- asum + apply(a[kyear,,],2,sum,na.rm=TRUE)*apply(rho.est[,1:(kyear-1)],1,prod)
           }
           N1.est <- matrix(ysum/asum,nrow=nT,ncol=nS,byrow=TRUE)
           # Concentrared likelihood formulation
           llik <- array(NA,dim=c(nY,nT,nS))
           llik[1,,] <- dpois(y[1,,],N1.est*a[1,,],log=TRUE)
           llik[2,,] <- dpois(y[2,,],N1.est*a[2,,]*matrix(rho.est[,1],nrow=nT,ncol=nS,byrow=TRUE),log=TRUE)
           for(kyear in 3:nY){
             llik[kyear,,] <- dpois(y[kyear,,],N1.est*a[kyear,,]*matrix(apply(rho.est[,1:(kyear-1)],1,prod),nrow=nT,ncol=nS,byrow=TRUE),log=TRUE)}},
         "2"={ # bivoltine species
           # productivity for each generation can be constant, time varying, or varying with one or more covariates (additively)
           switch(rho.m,
                  const={
                    rho.est <- list(matrix(exp(parm[par.index+1]),nrow=nS,ncol=nY),matrix(exp(parm[par.index+2]),nrow=nS,ncol=nY-1)); par.index <- par.index + 2},
                  indyear={
                    rho.est <- list(matrix(exp(parm[(par.index+1):(par.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE),matrix(exp(parm[(par.index+nY+1):(par.index+(nY-1)*2+1)]),nrow=nS,ncol=nY-1,byrow=TRUE)); par.index <- par.index + (nY-1)*2+1},
                  covs={
                    switch(rho1.m,
                           l1cov={
                             rho1.est <- exp(parm[par.index + 1] + parm[par.index + 2]*rho1_cov1); par.index <- par.index + 2},
                           l2cov={
                             rho1.est <- exp(parm[par.index + 1] + parm[par.index + 2]*rho1_cov1 + parm[par.index+3]*rho1_cov2); par.index <- par.index + 3})
                    switch(rho2.m,
                           l1cov={
                             rho2.est <- exp(parm[par.index + 1] + parm[par.index + 2]*rho2_cov1); par.index <- par.index + 2},
                           l2cov={
                             rho2.est <- exp(parm[par.index + 1] + parm[par.index + 2]*rho2_cov1 + parm[par.index+3]*rho2_cov2); par.index <- par.index + 3})
                    rho.est <- list(rho1.est,rho2.est)})
           
           # mu for each generation can be constant, time varying, or varying with one or more covariates (additively)
           switch(mu.m,
                  const={
                    mu.est <- array(exp(parm[(par.index+1):(par.index+2)]),c(2,nS,nY)); par.index <- par.index+2},
                  indyear={
                    mu.est <- array(NA,c(2,nS,nY))
                    mu.est[1,,] <- matrix(exp(parm[(par.index+1):(par.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE)
                    mu.est[2,,] <- matrix(exp(parm[(par.index+nY+1):(par.index+2*nY)]),nrow=nS,ncol=nY,byrow=TRUE); par.index <- par.index+2*nY},
                  covs={
                    switch(mu1.m,
                           l1cov={
                             mu1.est <- exp(parm[par.index+1] + parm[par.index + 2]*mu1_cov1); par.index <- par.index + 2},
                           l2cov={
                             mu1.est <- exp(parm[par.index + 1] + parm[par.index + 2]*mu1_cov1 + parm[par.index+3]*mu1_cov2); par.index <- par.index + 3})
                    switch(mu2.m,
                           l1cov={
                             mu2.est <- exp(parm[par.index + 1] + parm[par.index + 2]*mu2_cov1); par.index <- par.index + 2},
                           l2cov={
                             mu2.est <- exp(parm[par.index + 1] + parm[par.index + 2]*mu2_cov1 + parm[par.index+3]*mu2_cov2); par.index <- par.index + 3})	
                    mu.est <- array(NA,c(2,nS,nY))
                    mu.est[1,,] <- mu1.est
                    mu.est[2,,] <- mu2.est})					
           
           # sigma for each generation can be constant, time varying, or varying with one or more covariates (additively)
           switch(sigma.m,
                  const={
                    sigma.est <-  array(exp(parm[(par.index+1):(par.index+2)]),c(2,nS,nY)); par.index <- par.index+2},
                  indyear={
                    sigma.est <- array(NA,c(2,nS,nY))
                    sigma.est[1,,] <- matrix(exp(parm[(par.index+1):(par.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE)
                    sigma.est[2,,] <- matrix(exp(parm[(par.index+nY+1):(par.index+2*nY)]),nrow=nS,ncol=nY,byrow=TRUE); par.index <- par.index+2*nY},
                  covs={
                    switch(sigma1.m,
                           l1cov={
                             sigma1.est <- exp(parm[par.index+1] + parm[par.index + 2]*sigma1_cov1); par.index <- par.index + 2},
                           l2cov={
                             sigma1.est <- exp(parm[par.index + 1] + parm[par.index + 2]*sigma1_cov1 + parm[par.index+3]*sigma1_cov2); par.index <- par.index + 3})
                    switch(sigma2.m,
                           l1cov={
                             sigma2.est <- exp(parm[par.index + 1] + parm[par.index + 2]*sigma2_cov1); par.index <- par.index + 2},
                           l2cov={
                             sigma2.est <- exp(parm[par.index + 1] + parm[par.index + 2]*sigma2_cov1 + parm[par.index+3]*sigma2_cov2); par.index <- par.index + 3})	
                    sigma.est <- array(NA,c(2,nS,nY))
                    sigma.est[1,,] <- sigma1.est
                    sigma.est[2,,] <- sigma2.est})
           
           switch(atype,
                  "N"={
                    a <- array(NA,dim=c(nY,2,nT,nS))
                    for(kyear in 1:nY){
                      a[kyear,1,,] <- matrix(dnorm(rep(1:nT,each=nS),mu.est[1,,kyear],sigma.est[1,,kyear]),nrow=nT,ncol=nS,byrow=TRUE)
                      a[kyear,2,,] <- matrix(dnorm(rep(1:nT,each=nS),mu.est[1,,kyear]+mu.est[2,,kyear],sigma.est[2,,kyear]),nrow=nT,ncol=nS,byrow=TRUE)}},
                  "SO"={
                    # survival probability for each generation can be constant, time varying, or varying with one or more covariates (additively)
                    switch(phi.m,
                           const={
                             phi.est <- array(rep(expit(parm[(par.index+1):(par.index+2)]),each=nS*nY),c(nS,nY,2)); par.index <- par.index + 2},
                           indyear={
                             phi.est <- array(c(rep(expit(parm[(par.index+1):(par.index+nY)]),each=nS),rep(expit(parm[(par.index+nY+1):(par.index+2*nY)]),each=nS)),c(nS,nY,2)); par.index <- par.index+2*nY},
                           covs={
                             phi.est <- array(NA,c(nS,nY,2))
                             switch(phi1.m,
                                    l1cov={
                                      phi.est[,,1] <- expit(parm[par.index+1] + parm[par.index+2]*phi1_cov1); par.index <- par.index + 2},
                                    l2cov={
                                      phi.est[,,1] <- expit(parm[par.index+1] + parm[par.index+2]*phi1_cov1 + parm[par.index+3]*phi1_cov2); par.index <- par.index + 3})
                             switch(phi2.m,
                                    l1cov={
                                      phi.est[,,2] <- expit(parm[par.index+1] + parm[par.index+2]*phi2_cov1); par.index <- par.index + 2},
                                    l2cov={
                                      phi.est[,,2] <- expit(parm[par.index+1] + parm[par.index+2]*phi2_cov1 + parm[par.index+3]*phi2_cov2); par.index <- par.index + 3})})
                    
                    # Array for betta, describing the proportions of N arriving		
                    betta.est1 <- array(c(pnorm(1,mean=mu.est[1,,],sd=sigma.est[1,,]),pnorm(rep(2:(nT-1),each=nS*nY),mean=mu.est[1,,],sd=sigma.est[1,,])-pnorm(rep(1:(nT-2),each=nS*nY),mean=mu.est[1,,],sd=sigma.est[1,,]),1-pnorm(nT-1,mean=mu.est[1,,],sd=sigma.est[1,,])),c(nS,nY,nT))
                    betta.est2 <- array(c(pnorm(1,mean=mu.est[1,,]+mu.est[2,,],sd=sigma.est[2,,]),pnorm(rep(2:(nT-1),each=nS*nY),mean=mu.est[1,,]+mu.est[2,,],sd=sigma.est[2,,])-pnorm(rep(1:(nT-2),each=nS*nY),mean=mu.est[1,,]+mu.est[2,,],sd=sigma.est[2,,]),1-pnorm(nT-1,mean=mu.est[1,,]+mu.est[2,,],sd=sigma.est[2,,])),c(nS,nY,nT))
                    a1 <- betta.est1
                    a2 <- betta.est2
                    a <- array(NA,dim=c(2,nS,nY,nT))
                    a[1,,,] <- array(cbind(c(betta.est1[,,1]),sapply(2:nT,FUN=bettafunc,betta.est1,a1,phi.est[,,1])),c(nS,nY,nT))
                    a[2,,,] <- array(cbind(c(betta.est2[,,1]),sapply(2:nT,FUN=bettafunc,betta.est2,a2,phi.est[,,2])),c(nS,nY,nT))
                    a <- aperm(a,c(3,1,4,2))})
           
           a[,1,,][is.na(y)] <- NA	
           a[,2,,][is.na(y)] <- NA	
           
           # Estimate N_{i,1,1} (equation A2) for use in the concentrated likelihood
           ysum <- apply(y,3,sum,na.rm=TRUE)
           asum <- apply(a[1,1,,],2,sum,na.rm=TRUE) + rho.est[[1]][,1]*apply(a[1,2,,],2,sum,na.rm=TRUE) + (apply(a[2,1,,],2,sum,na.rm=TRUE)+apply(a[2,2,,],2,sum,na.rm=TRUE)*rho.est[[1]][,2])*rho.est[[1]][,1]*rho.est[[2]][,1]	
           for(kyear in 3:nY){
             asum <- asum + (apply(a[kyear,1,,],2,sum,na.rm=TRUE)+apply(a[kyear,2,,],2,sum,na.rm=TRUE)*rho.est[[1]][,kyear])*apply(rho.est[[1]][,1:(kyear-1)],1,prod)*apply(rho.est[[2]][,1:(kyear-1)],1,prod)}
           N1.est <- matrix(ysum/asum,nrow=nT,ncol=nS,byrow=TRUE)
           
           # Concentrated likelihood formulation
           llik <- array(NA,dim=c(nY,nT,nS))
           llik[1,,] <- dpois(y[1,,],N1.est*(a[1,1,,]+rho.est[[1]][,1]*a[1,2,,]),log=TRUE)
           llik[2,,] <- dpois(y[2,,],N1.est*(a[2,1,,]+a[2,2,,]*matrix(rho.est[[1]][,2],nrow=nT,ncol=nS,byrow=TRUE))*matrix(rho.est[[1]][,1]*rho.est[[2]][,1],nrow=nT,ncol=nS,byrow=TRUE),log=TRUE)
           for(kyear in 3:nY){
             llik[kyear,,] <- dpois(y[kyear,,],N1.est*(a[kyear,1,,]+a[kyear,2,,]*matrix(rho.est[[1]][,kyear],nrow=nT,ncol=nS,byrow=TRUE))*matrix(apply(rho.est[[1]][,1:(kyear-1)],1,prod)*apply(rho.est[[2]][,1:(kyear-1)],1,prod),nrow=nT,ncol=nS,byrow=TRUE),log=TRUE)}})	
  -1*sum(llik,na.rm=TRUE) 
}

# Starting values
start_val_func <- function(){
  switch(M,
         "1"={
           switch(rho.m,
                  const={rho.st <- log(sample(seq(.75,1.25,.05),1))},
                  indyear={rho.st <- log(rep(sample(seq(.75,1.25,.05),1),nY-1))},
                  l1cov={rho.st <- c(log(sample(seq(.75,1.25,.05),1)),0); if(is.null(rho1_cov1)) stop("rho1_cov1 not specified")},
                  l2cov={rho.st <- c(log(sample(seq(.75,1.25,.05),1)),0,0); if(is.null(rho1_cov1) | is.null(rho1_cov2)) stop("rho1_cov1 and/or rho1_cov2 not specified")},
                  stop("rho.m should be const, indyear, l1cov or l2cov when M=1"))
           
           switch(mu.m,
                  const={mu.st <- log(sample(1:nT,1))},
                  indyear={mu.st <- rep(log(sample(1:nT,1)),nY)},
                  l1cov={mu.st <- c(log(sample(1:nT,1)),0); if(is.null(mu1_cov1)) stop("mu1_cov1 not specified")},
                  l2cov={mu.st <- c(log(sample(1:nT,1)),0,0); if(is.null(mu1_cov1) | is.null(mu1_cov2)) stop("mu1_cov1 and/or mu1_cov2 not specified")},
                  stop("mu.m should be const, indyear, l1cov or l2cov when M=1"))
           
           switch(sigma.m, 
                  const={sigma.st <- log(sample(2:3,1))},
                  indyear={sigma.st <- rep(log(sample(2:3,1)),nY)},
                  l1cov={sigma.st <- c(log(sample(2:3,1)),0); if(is.null(sigma1_cov1)) stop("sigma1_cov1 not specified")},
                  l2cov={sigma.st <- c(log(sample(2:3,1)),0,0); if(is.null(sigma1_cov1) | is.null(sigma1_cov2)) stop("sigma1_cov1 and/or sigma1_cov2 not specified")},
                  stop("sigma.m should be const, indyear, l1cov or l2cov when M=1"))
           
           phi.st <- NULL	
           if(atype == "SO"){
             switch(phi.m,
                    const={phi.st <- logit(sample(seq(.3,.9,.1),1))},
                    indyear={phi.st <- logit(rep(sample(seq(.3,.9,.1),1),nY))},
                    l1cov={phi.st <- c(logit(sample(seq(.3,.9,.1),1)),0); if(is.null(phi1_cov1)) stop("phi1_cov1 not specified")},
                    l2cov={phi.st <- c(logit(sample(seq(.3,.9,.1),1)),0,0); if(is.null(phi1_cov1) | is.null(phi1_cov2)) stop("phi1_cov1 and/or phi1_cov2 not specified")},
                    stop("phi.m should be const, indyear, l1cov or l2cov when M=1"))}},
         "2"={
           switch(rho.m,
                  const={rho.st <- log(sample(seq(.75,1.25,.05),2))},
                  indyear={rho.st <- log(rep(sample(seq(.75,1.25,.05),1),2*(nY-1)+1))},
                  covs={
                    switch(rho1.m,
                           l1cov={rho.st <- c(log(sample(seq(.75,1.25,.05),1)),0); if(is.null(rho1_cov1)) stop("rho1_cov1 not specified")},
                           l2cov={rho.st <- c(log(sample(seq(.75,1.25,.05),1)),0,0); if(is.null(rho1_cov1) | is.null(rho1_cov2)) stop("rho1_cov1 and/or rho1_cov2 not specified")},
                           stop("rho1.m should be l1cov or l2cov when M=2 and rho.m=covs"))
                    switch(rho2.m,
                           l1cov={rho.st <- c(rho.st,log(sample(seq(.75,1.25,.05),1)),0); if(is.null(rho2_cov1)) stop("rho2_cov1 not specified")},
                           l2cov={rho.st <- c(rho.st,log(sample(seq(.75,1.25,.05),1)),0,0); if(is.null(rho2_cov1) | is.null(rho2_cov2)) stop("rho2_cov1 and/or rho2_cov2 not specified")},
                           stop("rho2.m should be l1cov or l2cov when M=2 and rho.m=covs"))},
                  stop("rho.m should be const, indyear, or covs when M=2"))	
           
           switch(mu.m,
                  const={mu.st <- log(sample(1:round(nT/2),2))},
                  indyear={mu.st <- log(rep(sample(1:round(nT/2),2),nY))},
                  covs={
                    switch(mu1.m,
                           l1cov={mu.st <- c(log(sample(1:round(nT/2),1)),0); if(is.null(mu1_cov1)) stop("mu1_cov1 not specified")},
                           l2cov={mu.st <- c(log(sample(1:round(nT/2),1)),0,0); if(is.null(mu1_cov1) | is.null(mu1_cov2)) stop("mu1_cov1 and/or mu1_cov2 not specified")},
                           stop("mu1.m should be l1cov or l2cov when M=2 and mu.m=covs"))
                    switch(mu2.m,
                           l1cov={mu.st <- c(mu.st,log(sample(1:round(nT/2),1)),0); if(is.null(mu2_cov1)) stop("mu2_cov1 not specified")},
                           l2cov={mu.st <- c(mu.st,log(sample(1:round(nT/2),1)),0,0); if(is.null(mu2_cov1) | is.null(mu2_cov2)) stop("mu2_cov1 and/or mu2_cov2 not specified")},
                           stop("mu2.m should be l1cov or l2cov when M=2 and mu.m=covs"))},
                  stop("mu.m should be const, indyear, or covs when M=2"))	
           
           switch(sigma.m, 
                  const={
                    sigma.st <- log(sample(2:3,2))},
                  indyear={
                    sigma.st <- log(rep(sample(2:3,2),nY))},
                  covs={
                    switch(sigma1.m,
                           l1cov={sigma.st <- c(log(sample(2:3,1)),0); if(is.null(sigma_cov1)) stop("sigma1_cov1 not specified")},
                           l2cov={sigma.st <- c(log(sample(2:3,1)),0,0); if(is.null(sigma1_cov1) | is.null(sigma1_cov2)) stop("sigma1_cov1 and/or sigma1_cov2 not specified")},
                           stop("sigma1.m should be l1cov or l2cov when M=2 and sigma.m=covs"))
                    switch(sigma2.m,
                           l1cov={sigma.st <- c(sigma.st,log(sample(2:3,1)),0); if(is.null(sigma2_cov1)) stop("sigma2_cov1 not specified")},
                           l2cov={sigma.st <- c(sigma.st,log(sample(2:3,1)),0,0); if(is.null(sigma2_cov1) | is.null(sigma2_cov2)) stop("sigma2_cov1 and/or sigma2_cov2 not specified")},
                           stop("sigma2.m should be l1cov or l2cov when M=2 and sigma.m=covs"))},
                  stop("sigma.m should be const, indyear, or covs when M=2"))
           
           phi.st <- NULL		
           if(atype == "SO"){
             switch(phi.m,
                    const={phi.st <- logit(rep(sample(seq(.3,.9,.1),1),2))},
                    indyear={phi.st <- logit(rep(sample(seq(.3,.9,.1),1),2*nY))},
                    covs={
                      switch(phi1.m,
                             l1cov={phi.st <- c(logit(sample(seq(.3,.9,.1),1)),0); if(is.null(phi1_cov1)) stop("phi1_cov1 not specified")},
                             l2cov={phi.st <- c(logit(sample(seq(.3,.9,.1),1)),0,0); if(is.null(phi1_cov1) | is.null(phi1_cov2)) stop("phi1_cov1 and/or phi1_cov2 not specified")},
                             stop("phi1.m should be l1cov or l2cov when M=2 and phi.m=covs"))
                      switch(phi2.m,
                             l1cov={phi.st <- c(phi.st,logit(sample(seq(.3,.9,.1),1)),0); if(is.null(phi2_cov1)) stop("phi2_cov1 not specified")},
                             l2cov={phi.st <- c(phi.st,logit(sample(seq(.3,.9,.1),1)),0,0); if(is.null(phi2_cov1) | is.null(phi2_cov2)) stop("phi2_cov1 and/or phi2_cov2 not specified")},
                             stop("phi2.m should be l1cov or l2cov when M=2 and phi.m=covs"))},
                    stop("phi.m should be const, indyear, or covs when M=2"))}})
  parm <- c(rho.st,mu.st,sigma.st,phi.st)
}	

# Wrapper for fitting dynamic models for multiple starts	
fit_dynamic <- function(){
  if(!(atype %in% c("N","SO"))) stop("Error: atype must be N or SO")
  if(!(M %in% 1:2)) stop("Error: M must be 1 or 2")
  
  temp.fit <- list(); temp.ll <- rep(NA,nstart)
  for(kstart in 1:nstart){
    temp.fit[[kstart]] <- try(fit_func(),silent=FALSE)
    while(is.na(temp.fit[[kstart]][[1]])){
      temp.fit[[kstart]] <-  try(fit_func(),silent=FALSE)}
    temp.ll[kstart] <- temp.fit[[kstart]]$ll.val}
  output <- list(temp.fit[[ min(c(1:nstart)[temp.ll==max(temp.ll,na.rm=T)])]],temp.fit,temp.ll)
}

# Fit the dynamic models	
fit_func <- function(){
  parm <- start_val_func()
  st1 <- proc.time()
  this.fit <- try(optim(parm,dynamic_ll_func,hessian=TRUE,method="L-BFGS-B",control=list(maxit=1e4,trace=1)),silent=TRUE)
  et1 <- proc.time()
  
  if(is.list(this.fit) & class(try(solve(this.fit$hessian),silent=TRUE)) != "try-error"){	
    # Model output
    out.index <- 0		
    rho.int <- rho.slope <- mu.int <- mu.slope <- sigma.int <- sigma.slope <- phi.out <- phi.int <- phi.slope <- betta.out  <- NULL	
    
    switch(M,
           "1"={
             switch(rho.m,
                    const={
                      rho.out <- matrix(exp(this.fit$par[out.index+1]),nrow=nS,ncol=nY-1);out.index <- out.index + 1},
                    indyear={
                      rho.out <- matrix(exp(this.fit$par[(out.index+1):(out.index+nY-1)]),nrow=nS,ncol=nY-1,byrow=TRUE); out.index <- out.index + nY-1},
                    l1cov={
                      rho.out <- exp(this.fit$par[out.index+1] + this.fit$par[out.index+2]*rho1_cov1); rho.int <- this.fit$par[out.index+1]; rho.slope <- this.fit$par[out.index+2]; out.index <- out.index + 2},
                    l2cov={
                      rho.out <- exp(this.fit$par[out.index+1] + this.fit$par[out.index+2]*rho1_cov1 + this.fit$par[out.index+3]*rho1_cov2); rho.int <- this.fit$par[out.index+1]; rho.slope <- this.fit$par[(out.index+2):(out.index+3)]; out.index <- out.index + 3})
             
             switch(mu.m,
                    const={
                      mu.out <- matrix(exp(this.fit$par[out.index+1]),nrow=nS,ncol=nY);out.index <- out.index+1},
                    indyear={
                      mu.out <- matrix(exp(this.fit$par[(out.index+1):(out.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE); out.index <- out.index+nY},
                    l1cov={
                      mu.out <- exp(this.fit$par[out.index+1] + this.fit$par[out.index+2]*mu1_cov1); mu.int <- this.fit$par[out.index+1]; mu.slope <- this.fit$par[out.index+2]; out.index <- out.index+2},
                    l2cov={
                      mu.out <- exp(this.fit$par[out.index+1] + this.fit$par[out.index+2]*mu1_cov1 + this.fit$par[out.index+3]*mu1_cov2);mu.int <- this.fit$par[out.index+1]; mu.slope <- this.fit$par[(out.index+2):(out.index+3)]; out.index <- out.index+3})
             
             switch(sigma.m,
                    const={
                      sigma.out <- matrix(exp(this.fit$par[out.index+1]),nrow=nS,ncol=nY); out.index <- out.index+1},
                    indyear={
                      sigma.out <-  matrix(exp(this.fit$par[(out.index+1):(out.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE); out.index <- out.index+nY},
                    l1cov={
                      sigma.out <- exp(this.fit$par[out.index+1] +  this.fit$par[out.index+2]*sigma1_cov1); sigma.int <- this.fit$par[out.index+1]; sigma.slope <- this.fit$par[out.index+2]; out.index <- out.index+2},
                    l2cov={
                      sigma.out <- exp(this.fit$par[out.index+1] +  this.fit$par[out.index+2]*sigma1_cov1 + this.fit$par[out.index+3]*sigma1_cov2); sigma.int <- this.fit$par[out.index+1]; sigma.slope <- this.fit$par[(out.index+2):(out.index+3)]; out.index <- out.index+3})		
             
             switch(atype,
                    "N"={
                      a.out <- array(NA,dim=c(nY,nT,nS))
                      for(kyear in 1:nY){
                        a.out[kyear,,] <- matrix(dnorm(rep(1:nT,nS),rep(mu.out[,kyear],each=nT),sigma.out[kyear]),nrow=nT,ncol=nS)}},
                    "SO"={
                      switch(phi.m,
                             const={
                               phi.out <-  matrix(expit(this.fit$par[out.index+1]),nrow=nS,ncol=nY); out.index <- out.index + 1},
                             indyear={
                               phi.out <- matrix(expit(this.fit$par[(out.index+1):(out.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE); out.index <- out.index + nY},
                             l1cov={
                               phi.out <- expit(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*phi1_cov1); phi.int <- this.fit$par[out.index+1]; phi.slope <- this.fit$par[out.index+2]; out.index <- out.index + 2},
                             l2cov={
                               phi.out <- expit(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*phi1_cov1 + this.fit$par[out.index + 3]*phi1_cov2); phi.int <- this.fit$par[out.index+1]; phi.slope <- this.fit$par[(out.index+2):(out.index+3)]; out.index <- out.index + 3})
                      
                      betta.out <- array(c(pnorm(1,mean=mu.out,sd=sigma.out),pnorm(rep(2:(nT-1),each=nS*nY),mean=rep(mu.out,length(2:(nT-1))),sd=rep(sigma.out,length(2:(nT-1))))-pnorm(rep(1:(nT-2),each=nS*nY),mean=rep(mu.out,length(1:(nT-2))),sd=rep(sigma.out,length(1:(nT-2)))),1-pnorm(nT-1,mean=mu.out,sd=sigma.out)),c(nS,nY,nT))
                      a.out <- betta.out
                      a.out <- array(cbind(c(betta.out[,,1]),sapply(2:nT,FUN=	bettafunc,betta.out,a.out,phi.out)),c(nS,nY,nT))
                      a.out <- aperm(a.out,c(2,3,1))})	
             
             af.out <- a.out
             a.out[is.na(y)] <- NA	
             ysum <- apply(y,3,sum,na.rm=TRUE)
             asum <- apply(a.out[1,,],2,sum,na.rm=TRUE) + apply(a.out[2,,],2,sum,na.rm=TRUE)*rho.out[,1]
             for(kyear in 3:nY){
               asum <- asum + apply(a.out[kyear,,],2,sum,na.rm=TRUE)*apply(rho.out[,1:(kyear-1)],1,prod)}
             N1.out <- matrix(ysum/asum,nrow=nT,ncol=nS,byrow=TRUE)},
           "2"={	
             switch(rho.m,
                    const={
                      rho.out <- list(matrix(exp(this.fit$par[out.index+1]),nrow=nS,ncol=nY),matrix(exp(this.fit$par[out.index+2]),nrow=nS,ncol=nY-1)); out.index <- out.index + 2},
                    indyear={
                      rho.out <- list(matrix(exp(this.fit$par[(out.index+1):(out.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE),matrix(exp(this.fit$par[(out.index+nY+1):(out.index+(nY-1)*2+1)]),nrow=nS,ncol=nY-1,byrow=TRUE)); 		out.index <- out.index + (nY-1)*2+1},
                    covs={
                      rho.int <- this.fit$par[out.index + 1]
                      switch(rho1.m,
                             l1cov={
                               rho1.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*rho1_cov1); rho.slope <- this.fit$par[out.index + 2]; out.index <- out.index + 2},
                             l2cov={
                               rho1.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*rho1_cov1 + this.fit$par[out.index+3]*rho1_cov2); rho.slope <- this.fit$par[(out.index + 2):(out.index+3)];out.index <- out.index + 3})
                      rho.int <- c(rho.int,this.fit$par[out.index + 1]) 
                      switch(rho2.m,
                             l1cov={
                               rho2.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*rho2_cov1); rho.slope <- c(rho.slope,this.fit$par[out.index + 2]); out.index <- out.index + 2},
                             l2cov={
                               rho2.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*rho2_cov1 + this.fit$par[out.index+3]*rho2_cov2); rho.slope <- c(rho.slope,this.fit$par[(out.index + 2):(out.index+3)]); out.index <- out.index + 3})
                      rho.out <- list(rho1.out,rho2.out)})
             
             switch(mu.m,
                    const={
                      mu.out <- array(exp(this.fit$par[(out.index+1):(out.index+2)]),c(2,nS,nY)); out.index <- out.index+2},
                    indyear={
                      mu.out <- array(NA,c(2,nS,nY))
                      mu.out[1,,] <- matrix(exp(this.fit$par[(out.index+1):(out.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE)
                      mu.out[2,,] <- matrix(exp(this.fit$par[(out.index+nY+1):(out.index+2*nY)]),nrow=nS,ncol=nY,byrow=TRUE); out.index <- out.index+2*nY},
                    covs={
                      mu.int <- this.fit$par[out.index + 1]
                      switch(mu1.m,
                             l1cov={
                               mu1.out <- exp(this.fit$par[out.index+1] + this.fit$par[out.index + 2]*mu1_cov1); mu.slope <- this.fit$par[out.index + 2]; out.index <- out.index + 2},
                             l2cov={
                               mu1.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*mu1_cov1 + this.fit$par[out.index+3]*mu1_cov2); mu.slope <- this.fit$par[(out.index + 2):(out.index+3)]; out.index <- out.index + 3})
                      mu.int <- c(mu.int,this.fit$par[out.index + 1])
                      switch(mu2.m,
                             l1cov={
                               mu2.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*mu2_cov1); mu.slope <- c(mu.slope,this.fit$par[out.index + 2]); out.index <- out.index + 2},
                             l2cov={
                               mu2.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*mu2_cov1 + this.fit$par[out.index+3]*mu2_cov2); mu.slope <- c(mu.slope,this.fit$par[(out.index + 2):(out.index+3)]); out.index <- out.index + 3})	
                      mu.out <- array(NA,c(2,nS,nY))
                      mu.out[1,,] <- mu1.out
                      mu.out[2,,] <- mu2.out})					
             
             
             switch(sigma.m,
                    const={
                      sigma.out <-  array(exp(this.fit$par[(out.index+1):(out.index+2)]),c(2,nS,nY)); out.index <- out.index+2},
                    indyear={
                      sigma.out <- array(NA,c(2,nS,nY))
                      sigma.out[1,,] <- matrix(exp(this.fit$par[(out.index+1):(out.index+nY)]),nrow=nS,ncol=nY,byrow=TRUE)
                      sigma.out[2,,] <- matrix(exp(this.fit$par[(out.index+nY+1):(out.index+2*nY)]),nrow=nS,ncol=nY,byrow=TRUE); out.index <- out.index+2*nY},
                    covs={
                      sigma.int <- this.fit$par[out.index + 1]
                      switch(sigma1.m,
                             l1cov={
                               sigma1.out <- exp(this.fit$par[out.index+1] + this.fit$par[out.index + 2]*sigma1_cov1); sigma.slope <- this.fit$par[out.index + 2]; out.index <- out.index + 2},
                             l2cov={
                               sigma1.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*sigma1_cov1 + this.fit$par[out.index+3]*sigma1_cov2); sigma.slope <- this.fit$par[(out.index + 2):(out.index+3)]; out.index <- out.index + 3})
                      sigma.int <- c(sigma.int,this.fit$par[out.index + 1])
                      switch(sigma2.m,
                             l1cov={
                               sigma2.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*sigma2_cov1); sigma.slope <- c(sigma.slope,this.fit$par[out.index + 2]); out.index <- out.index + 2},
                             l2cov={
                               sigma2.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*sigma2_cov1 + this.fit$par[out.index+3]*sigma2_cov2); sigma.slope <- c(sigma.slope,this.fit$par[(out.index + 2):(out.index+3)]); out.index <- out.index + 3})	
                      sigma.out <- array(NA,c(2,nS,nY))
                      sigma.out[1,,] <- sigma1.out
                      sigma.out[2,,] <- sigma2.out})
             
             
             switch(atype,
                    "N"={
                      a.out <- array(NA,dim=c(nY,2,nT,nS))
                      for(kyear in 1:nY){
                        a.out[kyear,1,,] <- matrix(dnorm(rep(1:nT,each=nS),mu.out[1,,kyear],sigma.out[1,,kyear]),nrow=nT,ncol=nS,byrow=TRUE)
                        a.out[kyear,2,,] <- matrix(dnorm(rep(1:nT,each=nS),mu.out[1,,kyear]+mu.out[2,,kyear],sigma.out[2,,kyear]),nrow=nT,ncol=nS,byrow=TRUE)}},	
                    "SO"={
                      switch(phi.m,
                             const={
                               phi.out <- array(rep(expit(this.fit$par[(out.index+1):(out.index+2)]),each=nS*nY),c(nS,nY,2)); out.index <- out.index + 2},
                             indyear={
                               phi.out <- array(c(rep(expit(this.fit$par[(out.index+1):(out.index+nY)]),each=nS),rep(expit(this.fit$par[(out.index+nY+1):(out.index+2*nY)]),each=nS)),c(nS,nY,2)); out.index <- out.index+2*nY},
                             covs={
                               phi.out <- array(NA,c(nS,nY,2))
                               phi.int <- this.fit$par[out.index+1]
                               switch(phi1.m,
                                      l1cov={
                                        phi.out[,,1] <- expit(this.fit$par[out.index+1] + this.fit$par[out.index+2]*phi1_cov1); phi.slope <- this.fit$par[out.index + 2]; out.index <- out.index + 2},
                                      l2cov={
                                        phi.out[,,1] <- expit(this.fit$par[out.index+1] + this.fit$par[out.index+2]*phi1_cov1 + this.fit$par[out.index+3]*phi1_cov2); phi.slope <- this.fit$par[(out.index + 2):(out.index + 3)]; out.index <- out.index + 3})
                               phi.int <- c(phi.int,this.fit$par[out.index+1])
                               switch(phi2.m,
                                      l1cov={
                                        phi.out[,,2] <- expit(this.fit$par[out.index+1] + this.fit$par[out.index+2]*phi2_cov1); phi.slope <- c(phi.slope,this.fit$par[out.index + 2]); out.index <- out.index + 2},
                                      l2cov={
                                        phi.out[,,2] <- expit(this.fit$par[out.index+1] + this.fit$par[out.index+2]*phi2_cov1 + this.fit$par[out.index+3]*phi2_cov2); phi.slope <- c(phi.slope,this.fit$par[(out.index + 2):(out.index + 3)]); out.index <- out.index + 3})})
                      
                      a.out <- betta.out <- array(NA,dim=c(nY,2,nS,nT))
                      for(kyear in 1:nY){
                        betta.out[kyear,1,,] <- matrix(c(pnorm(1,mean=mu.out[1,,kyear],sd=sigma.out[1,,kyear]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.out[1,,kyear],length(2:(nT-1))),sd=rep(sigma.out[1,,kyear],length(2:(nT-1))))-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.out[1,,kyear],length(1:(nT-2))),sd=rep(sigma.out[1,,kyear],length(1:(nT-2)))),1-pnorm(nT-1,mean=mu.out[1,,kyear],sd=sigma.out[1,,kyear])),nrow=nS,ncol=nT)
                        betta.out[kyear,2,,] <- matrix(c(pnorm(1,mean=mu.out[1,,kyear]+mu.out[2,,kyear],sd=sigma.out[2,,kyear]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.out[1,,kyear]+mu.out[2,,kyear],length(2:(nT-1))),sd=rep(sigma.out[2,,kyear],length(2:(nT-1))))-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.out[1,,kyear]+mu.out[2,,kyear],length(1:(nT-2))),sd=rep(sigma.out[2,,kyear],length(1:(nT-2)))),1-pnorm(nT-1,mean=mu.out[1,,kyear]+mu.out[2,,kyear],sd=sigma.out[2,,kyear])),nrow=nS,ncol=nT)
                        
                        a.out[kyear,1,,] <- betta.out[kyear,1,,]
                        a.out[kyear,2,,] <- betta.out[kyear,2,,]
                        
                        for(j in 2:nT){
                          for(b in 1:(j-1)){			
                            a.out[kyear,1,,j] <- a.out[kyear,1,,j] + betta.out[kyear,1,,b]*phi.out[,kyear,1]^length(b:(j-1))
                            a.out[kyear,2,,j] <- a.out[kyear,2,,j] + betta.out[kyear,2,,b]*phi.out[,kyear,2]^length(b:(j-1))
                          }
                        }
                      }	
                      a.out <- aperm(a.out,c(1,2,4,3))})
             af.out <- a.out
             a.out[,1,,][is.na(y)] <- NA	
             a.out[,2,,][is.na(y)] <- NA	
             
             ysum <- apply(y,3,sum,na.rm=TRUE)
             asum <- apply(a.out[1,1,,],2,sum,na.rm=TRUE) + rho.out[[1]][,1]*apply(a.out[1,2,,],2,sum,na.rm=TRUE)	
             for(kyear in 2:nY){
               if(kyear==2){
                 asum <- asum + (apply(a.out[kyear,1,,],2,sum,na.rm=TRUE)+apply(a.out[kyear,2,,],2,sum,na.rm=TRUE)*rho.out[[1]][,kyear])*rho.out[[1]][,1]*rho.out[[2]][,1]
               } else {
                 asum <- asum + (apply(a.out[kyear,1,,],2,sum,na.rm=TRUE)+apply(a.out[kyear,2,,],2,sum,na.rm=TRUE)*rho.out[[1]][,kyear])*apply(rho.out[[1]][,1:(kyear-1)],1,prod)*apply(rho.out[[2]][,1:(kyear-1)],1,prod)}}
             N1.out <- matrix(ysum/asum,nrow=nT,ncol=nS,byrow=TRUE)})
    
    switch(M,
           "1"={
             N.out <- matrix(NA,nrow=nS,ncol=nY)
             N.out[,1] <- N1.out[1,]
             for(kyear in 2:nY){
               N.out[,kyear] <- N.out[,kyear-1]*rho.out[,kyear-1]
             }
             indexG <- apply(N.out,2,mean)
             FittedVal <-  array(apply(t(N.out),2,rep,26),dim(af.out))*af.out},
           "2"={
             N.out <- array(NA,c(2,nS,nY))
             N.out[1,,1] <- N1.out[1,]
             N.out[2,,1] <- N.out[1,,1]*rho.out[[1]][,1]
             for(kyear in 2:nY){
               N.out[1,,kyear] <- N.out[2,,kyear-1]*rho.out[[2]][,kyear-1]
               N.out[2,,kyear] <- N.out[1,,kyear-1]*rho.out[[1]][,kyear]}
             indexG <- matrix(c(apply(N.out[1,,],2,mean),apply(N.out[2,,],2,mean)),nrow=2,ncol=nY,byrow=TRUE)
             FittedVal <- array(apply(t(N.out[1,,]),2,rep,26),c(nY,nT,nS))*af.out[,1,,] + array(apply(t(N.out[2,,]),2,rep,26),c(nY,nT,nS))*af.out[,2,,]})
    
    
    Deviance <- 2*(sum((y*log(y/FittedVal)-(y-FittedVal))[!is.na(y) & y!=0])+sum(FittedVal[y== 0 & !is.na(y)]))
    ScaledDeviance <- Deviance/(length(y[!is.na(y)])-length(this.fit$par))
    
    
    
    output <- list(ll.val = -this.fit$value,
                   npar = length(this.fit$par),
                   modelspec = list(M = M,
                                    atype = atype,
                                    nS = nS,
                                    nT = nT,
                                    nY = nY,
                                    rho.m = rho.m,
                                    rho1.m = rho1.m,
                                    rho2.m = rho2.m,
                                    mu.m = mu.m,
                                    mu1.m = mu1.m,
                                    mu2.m = mu2.m,
                                    sigma.m = sigma.m,
                                    phi.m = phi.m,
                                    phi1.m = phi1.m,
                                    phi2.m = phi2.m),
                   modelresults = list(rho.out = rho.out,
                                       mu.out = mu.out,
                                       sigma.out = sigma.out,
                                       phi.out = phi.out,
                                       betta.out = betta.out,
                                       af.out = af.out,
                                       a.out = a.out,
                                       N1.out = N1.out,
                                       N.out = N.out,
                                       indexG = indexG,
                                       FittedVal = FittedVal,
                                       Deviance = Deviance,
                                       ScaledDeviance = ScaledDeviance,
                                       param = list(rho.int = rho.int, rho.slope = rho.slope, mu.int = mu.int, mu.slope = mu.slope, sigma.int = sigma.int, sigma.slope = sigma.slope, phi.int = phi.int,phi.slope = phi.slope)),
                   modelfit = list(convergence = this.fit$convergence,
                                   hessian = this.fit$hessian,
                                   allval = this.fit$par,			
                                   time = et1 - st1,
                                   starts = parm),
                   modeldata = list(y = y,
                                    covs = list(rho1_cov1 = rho1_cov1, rho2_cov1 = rho2_cov1, rho1_cov2 = rho1_cov2, rho2_cov2 = rho2_cov2, mu1_cov1 = mu1_cov1, mu1_cov2 = mu2_cov2, mu1_cov2 = mu1_cov2, mu2_cov2 = mu2_cov2, sigma1_cov1 = sigma1_cov1, sigma2_cov1 = sigma2_cov1, sigma1_cov2 = sigma1_cov2, sigma2_cov2 = sigma2_cov2, phi1_cov1 = phi1_cov1, phi2_cov1 = phi2_cov1, phi1_cov2 = phi1_cov2, phi2_cov2 = phi2_cov2)))
  } else {NA}
}
