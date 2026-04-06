mlag <- function(X,lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)  
}


hamiltonfilter <- function(beta1,beta2,sigma1,sigma2,p,q,y,x){
  T <- nrow(y)
  isig1 <- solve(sigma1)
  isig2 <- solve(sigma2)
  dsig1 <- sigma1^2
  dsig2 <- sigma2^2
  #Initialise the Filter
  pr_tr <- rbind(c(p[1],1-q[1]),c(1-p[1],q[1]))
  A <- rbind(diag(2)-pr_tr,matrix(1,1,2))
  EN <- matrix(c(0,0,1),3,1)
  ett <- ginv(crossprod(A))%*%t(A)%*%EN
  if (any(is.na(ett))){
    ett <- c(0.5,0.5)
  }
  #Forward Filter
  lik <- 0
  fprob <- matrix(0,T,2)
  for (it in 1:T){
    pr_tr <- rbind(c(p[it],1-q[it]),c(1-p[it],q[it]))
    em1 <- (y[it,] - x[it,]%*%beta1)
    em2 <- (y[it,] - x[it,]%*%beta2)
    
    neta1 <- (1/sqrt(dsig1))*exp(-0.5*(em1%*%isig1%*%t(em1))) 
    neta2 <- (1/sqrt(dsig2))*exp(-0.5*(em2%*%isig2%*%t(em2))) 
    
    neta1[neta1 < 1e-30] <- 1e-30
    neta2[neta2 < 1e-30] <- 1e-30
    neta1[neta1 > 1e30] <- 1e30
    neta2[neta2 > 1e30] <- 1e30
    
    #  if ((neta1 && neta2)==0) neta1 <- 0.1
    
    #Kim and Nelson Algorithm
    #  if (any(is.na(ett))) ett <- c(0.5,0.5)
    
    ett1 <- ett*rbind(neta1,neta2)
    
    fit <- sum(ett1)
    ett <- (pr_tr%*%ett1)/fit
    
    fprob[it,] <- t(ett1/fit)
    
    if (fit>0 && !is.na(fit)){
      lik <- lik+log(fit)
    }else{
      lik <- lik-10
    }
  }
  
  return(list(fprob=fprob,lik=lik))
}

getS <- function(fprob,p,q,Ncrit,maxdraws){
  PROBLEM <- 0
  j <- 1
  chck <- -1
  while(j<maxdraws && chck<0){
    ST1 <- getst(fprob,p,q)
    ST <- ST1[[1]]
    
    check1 <- length(ST[!duplicated(ST)])==2
    T1 <- sum(ST==0)
    T2 <- sum(ST==1)
    check2 <- T1>=Ncrit
    check3 <-  T2>=Ncrit# T2 <30#CHGGGGGGGG
    if (check1+check2+check3==3){
      chck=10
    }else{
      j=j+1
    }
  }
  if (check1+check2+check3<3){
    PROBLEM=1
  }
  return(list(ST=ST,prob=PROBLEM,smooth=ST1[[2]]))
}

getst <- function(fprob,p,q){
  T <- nrow(fprob)
  ST <- matrix(0,T,1)
  filtered <- matrix(0,T,1)
  
  p00 <- fprob[T,1]
  p01 <- fprob[T,2]
  r <- runif(1,0,1)
  ST[T,1] <- (r>=(p00/(p00+p01)))
  #go backwards
  for (it in (T-1):1){
    pr_tr <- rbind(c(p[it],1-q[it]),c(1-p[it],q[it]))
    if (ST[it+1]==0){
      p00 <- pr_tr[1,1]*fprob[it,1]
      p01 <- pr_tr[1,2]*fprob[it,2]
    }else if (ST[it+1]==1){
      p00 <- pr_tr[2,1]*fprob[it,1]
      p01 <- pr_tr[2,2]*fprob[it,2]
    }
    #sample regime numbers
    r <- runif(1,0,1)
    alph <- (p00/(p00+p01))
    if (r<alph){
      ST[it] <- 0
    }else{
      ST[it] <- 1
    }
    filtered[it] <- alph
  }
  
  return(list(ST,filtered))
}


KF_R <- function(y, Z,Ht,Qtt, m, p, t, B0, V0){
  #Define everything calculated down there
  # B0 <- as.vector(B0)
  # V0 <- as.matrix(V0)
  # Qtt <- as.matrix(Qtt)
  # bt <- matrix(0, t,m)
  # Vt <- matrix(0, m^2, t)
  # R <- matrix(0, p,m)
  # H <- matrix(0, t*m, p)
  # cfe <- matrix(0, p,1)
  # yt <- matrix(0, p, 1)
  # f <- matrix(0, p,p)
  # inv_f <- matrix(0, p,p)
  # btt <- matrix(0, m, 1)
  # Vtt <- matrix(0, m, m)
  
  bp <- B0 #the prediction at time t=0 is the initial state
  Vp <- V0 #Same for the variance
  bt <- matrix(0,t,m) #Create vector that stores bt conditional on information up to time t
  Vt <- matrix(0,m^2,t) #Same for variances
  
  for (i in 1:t){
    R <- Ht[i,] #CHK LATER NOISE INNOVATION VARIANCE TO MEASUREMENT ERRORS
    if(ncol(Qtt) > 1) Qt <- diag(Qtt[i,]) else Qt <- Qtt[i,]
    H <- Z#[i,,drop = F]
    
    cfe <- y[i] - H%*%bp   # conditional forecast error
    f <- H%*%Vp%*%t(H) + R    # variance of the conditional forecast error
    inv_f <- try(t(H)%*%solve(f), silent = T)
    if(is(inv_f, "try-error")) inv_f <- t(H)%*%ginv(f)
    btt <- bp + Vp%*%inv_f%*%cfe  #updated mean estimate for btt Vp * inv_F is the Kalman gain
    Vtt <- Vp - Vp%*%inv_f%*%H%*%Vp #updated variance estimate for btt
    if (i < t){
      bp <- btt
      Vp <- Vtt + Qt
    }
    bt[i,] <- t(btt)
    Vt[,i] <- matrix(Vtt,m^2,1)
  }
  
  # draw the final value of the latent states using the moments obtained from the KF filters' terminal state
  bdraw <- matrix(0,t,m)
  
  bdraw.temp <- try(btt+t(chol(Vtt))%*%rnorm(nrow(Vtt)), silent=T)
  if (is(bdraw.temp, "try-error")) bdraw.temp <- mvrnorm(1, btt, Vtt+diag(1e-6,m))
  bdraw[t,] <- bdraw.temp
  
  #Now do the backward recurssions
  for (i in 1:(t-1)){
    if(ncol(Qtt) > 1) Qt <- diag(Qtt[t-1,]) else Qt <- Qtt[t-1,]
    bf <- t(bdraw[t-i+1,])
    btt <- t(bt[t-i,])
    Vtt <- matrix(Vt[,t-i,drop=FALSE],m,m)
    f <- Vtt + Qt
    
    inv_f <- try(Vtt%*%solve(f), silent = T)
    if(is(inv_f, "try-error")) inv_f <- Vtt%*%ginv(f)
    
    cfe <- bf - btt
    bmean <- t(btt) + inv_f%*%t(cfe)
    bvar <- Vtt - inv_f%*%Vtt
    
    bdraw.temp <- try(bmean+t(chol(bvar))%*%rnorm(nrow(bvar)), silent=T)
    if (is(bdraw.temp, "try-error")) bdraw.temp <- mvrnorm(1, bmean, bvar+diag(1e-6,m))
    
    bdraw[t-i,] <- bdraw.temp
  }
  
  return(bdraw)
}


get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}
