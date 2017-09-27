##################################################################
# This file contains code to simulate simple example data to 
# demonstrate fitting the dynamic models.
##################################################################
# We show how to fit the dynamic model to univoltine data with 
# and without covariates for productivity
# A similar approach is taken for bivoltine data and including 
# covariates in mu, sigma and phi
##################################################################

# Set working directory, containing dynamic_model.r 
rm(list = ls())
setwd("C:/Users/boa10gb/Documents/R/Climate")
# Load file containing functions to fit the model
source("dynamic_model.r")

############################################
# We simulate data for the N_1 dynamic model
############################################
# Number of sites
nS <- 2
# Number of visits
nT <- 120
# Number of years
nY <- 6

# Set parameter values (here we assume the seasonal pattern to be the same across different sites and years) 
mu.t <- 60
sigma.t <- 10
a.t <- dnorm(1:nT,mu.t,sigma.t)
# Initial abundance values
N1.t <- rpois(nS,150)
# Productivity e.g. here we assume it varies with year but not site
rho.t <- seq(.75,1.5,length.out=nY-1)

# Simulate data
y <- array(NA,c(nY,nT,nS))
y[1,,] <- matrix(rpois(nT*nS,rep(N1.t,each=nT)*a.t),nrow=nT,ncol=nS)
y[2,,] <- matrix(rpois(nT*nS,rep(N1.t,each=nT)*rho.t[1]*a.t),nrow=nT,ncol=nS)
for(kyear in 3:nY){
  y[kyear,,] <- matrix(rpois(nT*nS,rep(N1.t,each=nT)*prod(rho.t[1:(kyear-1)])*a.t),nrow=nT,ncol=nS)}
# Create missing values		
y[sample(1:(nY*nT*nS),.3*nY*nT*nS)] <- NA

################################################
# Model fitting
################################################
# Example 1: productivity varying with each year
################################################
rho1.m <- rho2.m <- phi.m <- phi1.m <- phi2.m <- mu1.m <- mu2.m <- rho1_cov1 <- rho2_cov1 <- rho1_cov2 <- rho2_cov2 <- mu1_cov1 <- mu2_cov1 <- mu1_cov2 <- mu2_cov2 <- sigma1_cov1 <- sigma2_cov1 <- sigma1_cov2 <- sigma2_cov2 <-  phi1_cov1 <- phi2_cov1 <- phi1_cov2 <- phi2_cov2 <- NULL
# atype can be "N" or "SO"
atype <- "N"
# Specify the number of generations per year (i.e. univoltine/bivoltine)
M <- 1
# Specify formulation for rho, mu and sigma, as well as phi is atype == "SO"
rho.m <- "indyear"
mu.m <- "const"
sigma.m <- "const"
# Specify number of random starts 
nstart <- 3
# Fit the specfied dynamic model
output <- fit_dynamic()
# The output consists of a list of length 3: the best model of multiple starts, output from the multiple starts (as a list), and the log-likelihood values for the multiple starts

###################################################
# Example 2: productivity varies with a covariate
###################################################
# When included, covariates must always be in a nS x nY matrix, with the exception of rho1 in the univoltine case and rho2 in the bivoltine case, which should be a nS x (nY-1) matrix
# Scaling covariates is suggested.

rho1.m <- rho2.m <- phi.m <- phi1.m <- phi2.m <- mu1.m <- mu2.m <- rho1_cov1 <- rho2_cov1 <- rho1_cov2 <- rho2_cov2 <- mu1_cov1 <- mu2_cov1 <- mu1_cov2 <- mu2_cov2 <- sigma1_cov1 <- sigma2_cov1 <- sigma1_cov2 <- sigma2_cov2 <-  phi1_cov1 <- phi2_cov1 <- phi1_cov2 <- phi2_cov2 <- NULL
# atype can be "N" or "SO"
atype <- "N"
# Specify the number of generations per year (i.e. univoltine/bivoltine)
M <- 1
# Here we specify rho as varying with a single covariate
rho.m <- "l1cov"
mu.m <- "const"
sigma.m <- "const"
# Specify number of random starts 
nstart <- 3
# Covariate for rho
rho1_cov1 <- matrix(scale(1:(nY-1)),nrow=nS,ncol=nY-1,byrow=TRUE)
# Fit the specfied dynamic model
output <- fit_dynamic()

