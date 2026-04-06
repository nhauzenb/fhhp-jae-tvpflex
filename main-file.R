################################################################################
###----- General Bayesian time-varying parameter vector autoregressions -----###
###----------------- for modeling government bond yields -------------------####
################################################################################
rm(list = ls())
w.dir <- ""

###--------------------------------------------------------------------------###
###--------------------------- Preliminaries --------------------------------###
###--------------------------------------------------------------------------###

library(snow)
library(snowfall)
library(RcppArmadillo)
library(Rcpp)
library(coda)
library(stochvol)
library(MASS)
library(zoo)
library(invgamma)
library(scales)
library(Matrix)
library(sparseMVN)

source(paste0(w.dir, "aux_funcs/tvp_flex_func.R"))
source(paste0(w.dir, "aux_funcs/aux_files.R"))
sourceCpp(paste0(w.dir, "aux_funcs/kf.cpp"))

###--------------------------------------------------------------------------###
###----------------------- Define model specification -----------------------###
###--------------------------------------------------------------------------###

# Data 
data.set <- "US_ylds"       # Yields
trans <- "I0"               # Stationarity 
info.set <- "S"             # Small Nelson Siegel factors
end.in <- 2019 +11/12
freq   <- 12
var.set <- c("L", "S", "C") # Three variable included in "S" info set
target.var <- "L" # Dependent variable for eq.-by-eq. estimation (choose either "L", "S", "C")
# Model 
cons <- TRUE   # Intercept
sv   <- TRUE   # SV
xi   <- 1e-4   # Relatively tight on state equation variance
p    <- 3      # No. of lags

# Modifier variant
si.exo <- TRUE # Exogenous effect modfiers
q <- 1         # No. of RW factors
ms.fac <- TRUE # MS factor

# MCMC
nburn <- 1000
nsave <- 1000
thin <- 1

###--------------------------------------------------------------------------###
###--------------------------- Data setup -----------------------------------###
###--------------------------------------------------------------------------###
load(paste0(w.dir, "data/US_ylds.rda"))
if(trans == "I0"){
  Xraw <- apply(ylds[,2:16], 2, diff)
  Xraw.lev <- ylds[,2:16]
}else{
  Xraw <- Xraw.lev <- ylds[,2:16]
}
  
Xraw.sd <- apply(Xraw, 2, function(x){sd(x)})
Xraw <- apply(Xraw, 2, function(x){x/sd(x)})
Xraw <- ts(Xraw, end = c(2019,12), frequency = freq)
  
Xraw.lev.sd <- apply(Xraw.lev, 2, function(x){sd(x)})
Xraw.lev <- apply(Xraw.lev, 2, function(x){x/sd(x)})
Xraw.lev <- ts(Xraw.lev, end = c(2019,12), frequency = freq)
 
###----------------------- Nelson Siegel factors ----------------------------###
maturities <- c(1:15)*12
lambda <- 0.0609
Lambda <- matrix(0,length(maturities),3) #three Nelson-Siegel factors
for (jj in 1:length(maturities)){
    l11 <- (1-exp(-maturities[[jj]]*lambda))/(maturities[[jj]] * lambda)
    Lambda[jj,] <- c(1, l11, l11-exp(-maturities[[jj]]*lambda))
}

NS.fac.lev <- NS.fac <- ts(matrix(NA, nrow(Xraw), 3), start = time(Xraw)[1], frequency = 12)
  
for(tt in 1:nrow(Xraw))
{
  NS.fac[tt,] <- t(solve(crossprod(Lambda))%*%crossprod(Lambda,Xraw[tt,]))
  NS.fac.lev[tt,] <- t(solve(crossprod(Lambda))%*%crossprod(Lambda,Xraw.lev[tt,]))
}
  
colnames(NS.fac) <- colnames(NS.fac.lev) <- c("L", "S", "C")
ts.plot(NS.fac, type = "l")
rownames(Lambda) <- colnames(Xraw)
  
if(info.set == "S") Xraw <- NS.fac
  
###--------------------- External factors driving TVPs ----------------------###
if(si.exo){
    #Xraw.exo <- exo.df[,c("SMB", "HML", "PS_LEVEL", "NFCI")]
    Xraw.exo <- exo.df[,c("NFCI", "RF", "REC")]
    Xraw.exo[,1] <- apply(Xraw.exo[,1,drop = F], 2, function(x){(x-mean(x))/sd(x)})
    Xraw.exo <- ts(Xraw.exo, end = c(2019,12), frequency = freq)
    Xraw.exo <- window(Xraw.exo, start = start(NS.fac), end = end(NS.fac))
    si <- Xraw.exo
    si <- apply(si, 2, function(x)(x - mean(x))/sd(x))
    sl.z <- 1:ncol(si)
    
}else{
  si <- c()  
  sl.z <- NULL
}
  
# MS factor driving TVPs
if(ms.fac){
  si <- cbind(si, si.ms = rep(1, nrow(NS.fac)))
  sl.ms <- ncol(si) #indicator for switching indicator
  si.ms00 <- sample(x = 0:1, size = T, c(0.1, 0.9), replace = T)
}else{
  sl.ms <- NULL
}

##Latent factor driving TVPs
if(q > 0){
  lat.fac <- TRUE
  si.fac <- matrix(1, nrow(NS.fac), q)
  si <- cbind(si, si.fac)
  sl.fac <- ncol(si) - (q-1):0  # indictors for latent factors
}else{
  lat.fac <- FALSE
  sl.fac <- NULL
}

if(all(!si.exo, !lat.fac, !ms.fac)) si <- matrix(rep(1,nrow(NS.fac)))

G <- ncol(si)
si <- ts(si, start = start(NS.fac), frequency = 12)

sl.ind <- list(sl.z = sl.z, sl.ms = sl.ms, sl.fac = sl.fac)

# Covariances
if(target.var != var.set[[1]]){
  contemp.var <- 1:(which(var.set == target.var) -1)
}else{
  contemp.var <- NULL
}

Xraw <- Xraw[,var.set]
Xraw <- window(Xraw, end = end.in)
si <- as.matrix(window(si, end = end.in))
Yraw <- NULL

###--------------- Create specification list -------------------------------###
M <- ncol(Xraw)
T <- nrow(Xraw)
K <- M*p
K.v <- length(contemp.var)

model.setup <- list()
model.setup$target.var <- target.var
model.setup$contemp.var <- contemp.var
model.setup$var.set <- var.set
model.setup$ms.fac <- ms.fac
model.setup$p <- p
model.setup$cons <- cons
model.setup$ms.fac <- ms.fac
model.setup$lat.fac <- lat.fac
model.setup$q <- q
model.setup$sl.ind <- sl.ind
model.setup$G <- G
model.setup$sv <- sv
model.setup$xi <- xi

Regr.obj <- tvp.flex.func(Yraw = Yraw, Xraw = Xraw, thin = thin, nsave = nsave, nburn = nburn, model.setup = model.setup)
  


list2env(x = Regr.obj, envir = .GlobalEnv)
  
k.i <- ncol(Xraw)
if(KKK >= 1){
    A0.it_store <- array((-1)*bt.store[,,(KK+1):(KK+KKK),drop =F], c(nsave,T,KKK))
    A0.i_store <- matrix((-1)*b.store[,(KK+1):(KK+KKK),drop =F],nsave, KKK)
    
    mu0.i_store <- (-1)*mu.store[,(KK+1):(KK+KKK),,drop =F]
    hslam0.i_store <- lambda.mu_store[,(KK+1):(KK+KKK),,drop =F]
    
    dimnames(A0.it_store) <- list(NULL, NULL, var.set[contemp.var])
    dimnames(mu0.i_store) <- dimnames(hslam0.i_store) <- list(NULL, var.set[contemp.var], 1:G)
    
    colnames(A0.i_store) <- var.set[contemp.var]
}else{
    A0.it_store <- A0.i_store <- mu0.i_store <- hslam0.i_store <- NULL
}
  
if(cons){
    A.it_store <- bt.store[,,c(1:KK,K)]
    A.i_store <- b.store[,c(1:KK,K)]
    
    mu.i_store <- mu.store[,c(1:KK,K),]
    hslam.i_store <- lambda.mu_store[,c(1:KK,K),]
    
}else{
    A.it_store <- bt.store[,,c(1:KK)]
    A.i_store <- b.store[,c(1:KK)]
    
    mu.i_store <- mu.store[,1:KK,]
    hslam.i_store <- lambda.mu_store[,1:K,]
    
} 
  
res.list <- list(A.i_store = A.i_store, A0.i_store = A0.i_store, mu0.i_store = mu0.i_store, A.it_store = A.it_store, A0.it_store = A0.it_store, eht.store = eht.store, fit.store = fit.store, si.store = si.store, fprob.store  = fprob.store)

folder <- paste0(data.set, "_", info.set)
dir.create(folder)
save(res.list, file = paste0(folder, "/", target.var, ".rda"))  
