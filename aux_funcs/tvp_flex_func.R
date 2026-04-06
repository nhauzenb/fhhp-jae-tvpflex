###--------------------------------------------------------------------------###
###----------------------------- MCMC sampler -------------------------------### 
###---------------------- for flexible TVP regression -----------------------### 
###--------------------------------------------------------------------------###
tvp.flex.func <- function(Yraw = y, Xraw = X, thin = thin, nsave = nsave, nburn = nburn, model.setup = model.setup){

###--------------------------------------------------------------------------###
###-------------------- Extract specification -------------------------------###
###--------------------------------------------------------------------------###
list2env(model.setup,globalenv())

ntot <- nsave*thin + nburn
save.set <- seq(thin,nsave*thin,thin) + nburn
save.len <- length(save.set)
save.ind <- 0

###--------------------------------------------------------------------------###
###------------------------- Regression setup -------------------------------###
###--------------------------------------------------------------------------###
if (is.null(Yraw)){
   #Do the VAR lag part and compute y based on p lags
   k.i <- ncol(Xraw)
   X <- mlag(Xraw, p)
   si <- si[(p+1-1):(nrow(Xraw)-1),,drop =F]
   Yraw <- Xraw[(p+1):nrow(Xraw),,drop=F]
   X <- X[(p+1):nrow(X),]
   y <- Yraw[,target.var, drop = F]
   X.contemp <- Yraw[,contemp.var, drop = F]
   KKK <- ncol(X.contemp) #covariances
   KK <- k.i*p #coefficients
   
   if(cons) X <- cbind(X,X.contemp, 1) else  X <- cbind(X, X.contemp)
   
}else{
   y <- Yraw
   X <- Xraw
   KK <- ncol(X)
   KKK <- NULL
}

K <- ncol(X)
T <- nrow(y)
KG <- K + K*G 

###--------------------------------------------------------------------------###
###------------------- Default priors and Starting values -------------------###
###--------------------------------------------------------------------------###
b0 <- matrix(0,KG, 1)
V0 <- 1e-4*diag(KG)
V0inv <- diag(1/diag(V0))
Xnorm <- matrix(rep(as.vector(X), G+1), T, KG)
XpX <- crossprod(X)

# HS time-invariant part 
lambda.cons <- 1
nu.cons <- 1
tau.cons <- 1
zeta.cons <- 1

# HS on loadings 
lambda.mu <- matrix(1,K, G)
nu.mu <-  matrix(1,K, G)
tau.mu <- rep(1, G)
zeta.mu <- rep(1, G)
psi.mu <- matrix(NA, K, G)

# Variance 
## Homoskedastic
t0 <- 0.01
S0 <- 0.01 
## SV
if(sv){
  ## Priors for SV
  sv_priors <- specify_priors(
    mu = sv_normal(mean = 0, sd = 100), # prior on unconditional mean in the state equation
    phi = sv_beta(shape1 = 5, shape2 = 1.5), #informative prior to push the model towards a random walk in the state equation (for comparability)
    sigma2 = sv_gamma(shape = 0.5, rate = 0.5), # Gamma prior on the state innovation variance
    nu = sv_infinity(),
    rho = sv_constant(0))
  
  # Initialization of SV processes
  temp_sv <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
}

eht <- rep(1,T)
ht <- log(eht)

# Const part and factor loadings
obs.draw <- matrix(0,KG,1)
b.draw <- matrix(0, K, 1)
mu.mat <- matrix(0, K, G)
mu.draw <- as.vector(mu.mat)
psi_draw <- matrix(0.1,KG,1)

# TVP part
XXinv <- solve(crossprod(X))
diagXXinv <- diag(XXinv)
ub.sc <- xi*diagXXinv*T/K^2
Vd    <- ub.sc
Vmat  <- diag(ub.sc)

# Initial state factors
fac.b0 <- rep(0,q) 
fac.V0 <- diag(1e-10, q) # Fix initial state of RW factor to zero (see Chan et al. (2020))
fit.fac <- matrix(0, T, G)
Sig2.t <- matrix(0, T, 1) 
for(qq in 1:G){
  for(tt in 1:T){
    fit.fac[tt,qq] <-  as.numeric(X[tt,]%*%mu.mat[,qq,drop = F]%*%t(si[tt,qq,drop = F]))
    Sig2.t[tt,] <- t(X[tt,])%*%Vmat%*%X[tt,] + eht[tt]
  }  
}  

fprob <- rep(1, T)

###--------------------------------------------------------------------------###
###-------------------------------- Storage ---------------------------------###
###--------------------------------------------------------------------------###
# Coefficients
fit.store <- matrix(NA, nsave, T)
b.store <- omega.store <- matrix(NA,nsave,K)
eigen.store <- matrix(NA,nsave,T)
si.store <- array(NA, c(nsave, T, G))
mu.store <- array(NA,c(nsave,K,G))
fprob.store <- array(NA, c(nsave, T,2))
bt.store <- mu_t.store <- array(NA,c(nsave,T,K))
eht.store <- matrix(NA,nsave,T)
# Shrinkage parameters
psi_store <- matrix(NA,nsave,KG)
tau.cons_store <- matrix(NA,nsave,1)
lambda.cons_store <- matrix(NA,nsave,K)
tau.mu_store <- matrix(NA,nsave,G)
lambda.mu_store <- array(NA,c(nsave,K,G))

start <- Sys.time()
pb <- txtProgressBar(min = 0, max = ntot, style = 3) #start progress bar
###--------------------------------------------------------------------------###
###---------------------------- START MCMC loop -----------------------------###
###--------------------------------------------------------------------------###
for(irep in 1:ntot)
{
###--------------------------------------------------------------------------###
###---------------- Step 1: Sample time-invariant coefficients --------------###
###--------------------------------------------------------------------------###
for(ii in 1:T) Xnorm[ii,] <- c(X[ii,], kronecker(si[ii,], X[ii,]))  
norm.vec <- as.numeric(1/sqrt(Sig2.t))
Xnorm <- Xnorm*norm.vec
ynorm <- y*norm.vec
V1 <- try(solve(crossprod(Xnorm) + V0inv), silent=TRUE)
if (is(V1,"try-error")) V1 <- ginv(crossprod(Xnorm) + V0inv)
b1 <- V1%*%(crossprod(Xnorm,ynorm) + V0inv%*%b0)

obs.draw <- try(b1 + t(chol(V1))%*%rnorm(KG), silent=T)
if (is(obs.draw, "try-error")) obs.draw <- mvrnorm(1, b1, V1)
  
b.draw  <- obs.draw[1:K]
mu.draw <- obs.draw[(K+1):KG]
mu.mat  <- matrix(mu.draw, K, G)

###--------------------------------------------------------------------------###
###---------------- Step 2: Sample shrinkage parameters ---------------------###
###--------------------------------------------------------------------------###
#Step 2.1: Sample the prior for constant part using HS prior
hs_draw <- get.hs(bdraw=as.numeric(b.draw),lambda.hs=lambda.cons,nu.hs=nu.cons,tau.hs=tau.cons,zeta.hs=zeta.cons)
    
lambda.cons <- hs_draw$lambda
nu.cons <- hs_draw$nu
tau.cons <- hs_draw$tau
zeta.cons <- hs_draw$zeta
    
psi.cons <- hs_draw$psi
psi.cons[psi.cons<1e-20] <- 1e-20
psi.cons[psi.cons>100] <- 100

#Step 2.2: Sample the prior scaling factors for the loadings
for (k in seq_len(G)) {
    hs_mu <- get.hs(bdraw=as.numeric(mu.mat[,k]),lambda.hs=lambda.mu[,k],nu.hs=nu.mu[,k],tau.hs=tau.mu[k],zeta.hs=zeta.mu[k])
      
    lambda.mu[,k] <- hs_mu$lambda
    nu.mu[,k] <- hs_mu$nu
    tau.mu[k] <- hs_mu$tau
    zeta.mu[k] <- hs_mu$zeta
    psi.mu[,k] <- hs_mu$psi
}
psi.mu[psi.mu <1e-20] <- 1e-20
psi.mu[psi.mu > 100] <- 100

# Update prior variances
V0 <- diag(c(psi.cons, as.vector(psi.mu)))
V0inv <- diag(1/c(psi.cons, as.vector(psi.mu)))  

###--------------------------------------------------------------------------###
###--------------- Step 3: Sample latent states and TVPs --------------------###
###--------------------------------------------------------------------------###
Sig2.t <- matrix(0, T, 1) # joint variance
# Compute marginal fit and variance by substituting the state equation into the observation equation
for(qq in 1:G){
  for(tt in 1:T){
    fit.fac[tt,qq] <-  as.numeric(X[tt,]%*%mu.mat[,qq,drop = F]%*%t(si[tt,qq,drop = F]))
    Sig2.t[tt,]    <-  t(X[tt,])%*%Vmat%*%X[tt,] + eht[tt]
  }  
}  

# 3.1 Sample RW factors
if(lat.fac){ 
  y.fac.i <- y - X%*%b.draw - rowSums(fit.fac[,-sl.fac,drop = F]) 
  mu.fac.i <- X%*%mu.mat[,sl.fac, drop = F]
  si.fac.i <- si[,sl.fac,drop = F]
  Qt.i <- rep(1, q) # Fix state variance of factor to one (see Chan et al. (2020))
  
  # Draw dynamic factor with FFBS
  si.fac.i <- t(KF_fast(y=t(y.fac.i),Z=mu.fac.i,Ht=Sig2.t,Qtt=t(matrix(Qt.i, q, T)),m=q,p=1,t=T,B0=fac.b0,V0=fac.V0))
  # si.fac.i <- (si.fac.i - mean(si.fac.i))/sd(si.fac.i)
  si[,sl.fac] <- si.fac.i
}

# 3.2 Sample a single indicator using Markov switching process
if(ms.fac){
  y.ms <- (y - X%*%b.draw - rowSums(fit.fac[,-sl.ms, drop = F]))/sqrt(as.numeric(Sig2.t))
  si.ms <- si[,sl.ms]
  x.ms <-  X/sqrt(as.numeric(Sig2.t))
  
  trans <- table(factor(paste0(head(si.ms,-1),tail(si.ms,-1)), levels = c("00", "01", "10", "11")))
  p.draw <- rbeta(1, 10 + trans["11"], 1 + trans["10"])
  q.draw <- rbeta(1, 10 + trans["00"], 1 + trans["01"])
  ppnew <- rep(p.draw, T)
  qqnew <- rep(q.draw, T)
  
  mu.ms.1 <- matrix(rep(0,K))
  mu.ms.2 <- mu.mat[,sl.ms,drop = F]
  st_filt <- hamiltonfilter(beta1 = mu.ms.1, beta2 = mu.ms.2, sigma1 = 1, sigma2 = 1,ppnew,qqnew,y.ms,x.ms)
  fprob <- st_filt$fprob;liki <- st_filt$lik
  s_sample <- getS(fprob,ppnew,qqnew,10,20)
  problems <- s_sample$prob
  si.ms <- as.numeric(s_sample$ST)
  if (problems==1) si.ms <- si.ms00 else si.ms00 <- si.ms

  si[,sl.ms] <- si.ms
}

# Step 3.3: Define the TVPs conditional on the latent factors and everything else 
b.full <- state.mean <- matrix(NA, T, K)
cond.mean <- si%*%t(mu.mat)

for (tt in seq_len(T)){
  Qn <- 1/Sig2.t[tt] 
  eps.tt <- as.numeric(y[tt] - X[tt,]%*%b.draw - X[tt, ]%*%cond.mean[tt,])
  eta.mean <- Vd*X[tt,]*eps.tt*Qn 
  state.mean[tt,]  <- t(cond.mean[tt, ] + eta.mean) 
  
  VX <- Vd*X[tt,,drop=F]
  Var.state <- Vmat - t(VX)%*%Qn%*%VX
  
  Var.state.diag <- diag(Var.state)
  Var.state[abs(Var.state) < 1e-6] <- 0 #Because off-diagonals are extremely close to zero --> zero it out to speed up sampling
  diag(Var.state) <- Var.state.diag
  
  if (isDiagonal(Var.state)){
    b.full[tt,] <- b.draw + state.mean[tt,]+ rnorm(K, 0, sqrt(Var.state.diag))
  }else{
    Var.chol <- try(t(chol(Var.state)), silent = TRUE)
    if (is(Var.chol, "try-error")) b.full[tt,] <- b.draw + mvrnorm(1, state.mean[tt,], Var.state) else b.full[tt,] <- b.draw + state.mean[tt,] + Var.chol%*%rnorm(K)
  }
}

#Step 3.4: Sample scaling elements for the process innovation variances
shocks.states <- b.full - state.mean
SSE.states <- apply(shocks.states^2,2,sum)
for (ll in seq_len(K)) Vd[ll] <- GIGrvg::rgig(1, 1/2 - T/2, SSE.states[[ll]], 1/(2*1e-10))
Vd[Vd > ub.sc] <- ub.sc[Vd > ub.sc]
Vmat <- diag(Vd)

###--------------------------------------------------------------------------###
###--------------------- Step 4: Sample variances ---------------------------###
###--------------------------------------------------------------------------###
eta <- fit <- c()
for(tt in 1:T){
  fit[tt] <- X[tt,]%*%b.full[tt,]
  eta[tt] <- y[tt] - X[tt,]%*%b.full[tt,]
}
  
if (sv){ 
# Draw stochastic volatility parameters
  temp_sv <- svsample_fast_cpp(as.numeric(eta), startpara = temp_sv, startlatent = ht, priorspec = sv_priors)
  temp_sv[c("mu", "phi", "sigma", "nu", "rho")] <- as.list(temp_sv$para[, c("mu", "phi", "sigma", "nu", "rho")])
  ht <- t(temp_sv$latent)
  ht[ht < -12] <- -12 #off-set in case the model heavily overfits
  eht <- exp(as.numeric(ht))

}else{
    t1 <- (t0+T)/2
    S1 <- (S0 + sum(eta^2))/2
  
    sigma2 <-  as.numeric(1/rgamma(1,t1,as.numeric(S1)))
    eht   <- as.numeric(rep(sigma2,T))
}


###--------------------------------------------------------------------------###
###------------------------ Final step: Storage  ----------------------------###  
###--------------------------------------------------------------------------###
if(irep %in% save.set)
{
    save.ind <- save.ind + 1
    
    fit.store[save.ind,]  <- fit 
    
    b.store[save.ind,] <- as.numeric(b.draw)
    bt.store[save.ind,,] <- b.full
    eht.store[save.ind,] <- eht

    fprob.store[save.ind,,] <- fprob
    
    si.store[save.ind,,] <- si
    mu.store[save.ind,,] <- mu.mat
    mu_t.store[save.ind,,] <- si%*%t(mu.mat)
    
    omega.store[save.ind,] <- as.numeric(Vd)
    tau.cons_store[save.ind,] <- tau.cons
    lambda.cons_store[save.ind,] <- lambda.cons
    tau.mu_store[save.ind,] <- as.numeric(tau.mu)
    lambda.mu_store[save.ind,,] <- lambda.mu
    
}
setTxtProgressBar(pb, irep)
}

end <- Sys.time()
# Time in minutes
time.min <- (ts(end)-ts(start))/60


return(list(b.store = b.store, omega.store = omega.store, bt.store = bt.store, eht.store = eht.store, fit.store = fit.store, mu.store = mu.store, mu_t.store = mu_t.store, fprob.store = fprob.store, si.store = si.store, tau.cons_store = tau.cons_store, lambda.cons_store = lambda.cons_store, tau.mu_store = tau.mu_store, lambda.mu_store = lambda.mu_store, M = M, KK = KK, KKK= KKK, K = K, T = T, X = X, y =y,Xraw = Xraw, time.min = time.min, ntot = ntot))

}
