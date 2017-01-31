if("moments" %in% rownames(installed.packages()) == FALSE)
  install.packages("moments")
if("MASS" %in% rownames(installed.packages()) == FALSE)
  install.packages("MASS")
library("moments")
library("MASS")

## Psi Calculation ##
psi_cal <- function(x,est,p,q){
  mxpq <- max(p,q)
  nx <- length(x)
  psih <- rep(1,nx)
  pdim <- p+q+1
  omega <- as.matrix(est[1],ncol=1,nrow=1)
  alpha <- as.matrix(est[2:(1+p)],ncol=1,nrow=p)
  if(q==0){
    beta <- 0
  }else
  {
    beta <- as.matrix(est[(p+2):pdim],ncol=1,nrow=q) 
  }
  
  for( t in 1:mxpq){
    psih[t] <- est[1]
  }
  for( t in (mxpq+1):nx){
    pos1 <- t-1
    pos2 <- t-p
    pos3 <- t-q
    psih[t]=psit(x, psih, omega, alpha, beta, pos1, pos2, pos3, q)
  }
  return(psih)
}

## First Four Moments Calculation ##
eps_dist2 <- function(eps){
  
  mts <- all.moments(eps, order.max=4, central = TRUE)
  mue <- mean(eps)
  vare <- var(eps)
  skewe <- mts[4]/vare^(3/2)
  kurte <- mts[5]/vare^2
  dist1 <- "empirical"
  
  return(list(mue = mue, vare=vare, skewe=skewe, kurte=kurte, dist1=dist1))
}

## std error calculation ###
StdErrCal <- function(x,p,q,initval,momente){
  
  n <- length(x)
  pdim <- 1+p+q   
  info <-matrix(0,pdim,pdim)
  mue <- momente[1]
  vare <- momente[2]
  skewe <- momente[3]
  kurte <- momente[4]
  
  ## psi hat
  psih <- rep(1,n)
  
  mxpq <- max(p,q)
  for (t in 1:mxpq){
    psih[t] = initval[1]    
  }
  
  ## t=(mxpq+1):n  				
  for (t in (mxpq+1):n)
  {	
    omega <- as.matrix(initval[1],ncol=1,nrow=1)
    alpha <- as.matrix(initval[2:(1+p)],ncol=1,nrow=p)
    
    if(q != 0){
      beta <- as.matrix(initval[(p+2):pdim],ncol=1,nrow=q)
    }
    pos1 <- t-1
    pos2 <- t-p
    pos3 <- t-q
    
    psih[t]=psit(x, psih, omega, alpha, beta, pos1, pos2, pos3, q)
    if(abs(psih[t]) > log(.Machine$double.xmax/max(momente))/4) {
      psih[t] <- sign(psih[t])*log(.Machine$double.xmax/max(momente))/8
    }
    
    ## Define derivatives of psih(t) wrt theta: pdim*1 vector							
    ## First derivative									
    derpsi=derpsit(x, psih, pos1, pos2, pos3, pdim, q)
    
    ## Second derivative: pdim*pdim matrix									
    der2psi=der2psit(pdim)   
    
    ## mu(t), sigsq(t), gamma(t), kappa(t)							
    mu=mue*exp(psih[t]) 							
    sigsq=vare*exp(2*psih[t])						
    gamma=skewe*exp(3*psih[t]) 					
    kappa=kurte*exp(4*psih[t]) 
    
    ## First Derivatives of mu(t),sigsq(t), gamma(t) and kappa(t)												
    dermu=mue*exp(psih[t])*derpsi							
    dersigsq=2*vare*exp(2*psih[t])*derpsi
    
    ## Compute m(t) and M(t)							
    m = x[t]-mu       							
    qm = m**2-sigsq
    
    ## Compute Quadratic variations of m(t) and M(t) and 
    ## covariance of (m(t), M(t))							
    vm = sigsq
    vqm = kappa-vm**2  #Note: kappa is NOT excess kurtosis
    vmqm = gamma
    
    ## Define rho^2(t) 															
    termr = 1-((vmqm**2)/(vm*vqm))								
    rho = 1/termr
    if(is.na(rho) | is.infinite(rho))
      rho <- 1
    
    ## Define eta 								
    eta = vmqm/(vm*vqm)	
    if(is.na(eta) | is.infinite((eta)))
      eta = 0
    
    term1 <- (dermu%*%t(dermu))*(1/vm)
    term2 <- (dersigsq%*%t(dersigsq))*(1/vqm)
    term3 <- (dermu%*%t(dersigsq)+dersigsq%*%t(dermu))*eta
    info <- info + rho*(term1 + term2 - term3)
    
  }
  if(sum(is.na(info)) > 0 | sum(is.infinite(info)) > 0){
    stderr <- rep(1, pdim)
  }else{
    s <- svd(info)
    inf.index <- which(is.infinite(s$d))
    D <- diag(1/s$d)
    for(i in 1:length(s$d)){
      if(is.infinite(D[i,i])) 
        D[i,i] <- .Machine$double.xmin*100
    }
    info_inv <- s$v%*%D%*%t(s$u)
    
    stderr <- sqrt(abs(diag(info_inv)))
  }
  
  return(stderr)
}
########################################################################
### SIMULATION FUNCTION
########################################################################

## Simulate durations from the Log ACD1(p,q) model
## n: final sample size 
## nb: burn-in  
## Model parameters: omega, alpha, beta
## omega: scalar
## alpha = (alpha_1,alpha_2,...,alpha_p)
## beta = (beta_1,beta_2,...,beta_q)
## feps: distribution of the errors epsilon (positive valued), and can be 
##       exponential, weibull, or gamma
## par1 and par2: parameters of feps 
## default values for par1 and par2 (scale) are 1
## When feps is exponential, par1=lambda; we do not need par2
## When feps is weibull, par1=alpha and par2=beta =1
## When feps is gamma, par1=k (shape) and par2=theta (scale)=1
########################################################################


## version 2: based on Bowen and Giot ##
fsim.logacd <- function (n,nb,omega,alpha,beta,feps,par1=1,par2=1) 
{
  
  ########################################################################
  ## Check for stationarity constraint
  ########################################################################
  
  if (sum(alpha)+sum(beta)>=1) 
  { cat("Error: The simulated process is not weakly stationary")
  } else {
    
    ## Start simulating Durations
    nt <- n + nb  # Total number of simulated datapoints
    p <- length(alpha)
    q <- length(beta)
    mxpq <- max(p,q)		
    
    
    ########################################################################
    ## Initialize psis, xs, x, and psi vectors
    ########################################################################
    psis <- rep(1,nt)
    logxs <- rep(1,nt)
    xs <- rep(1,nt)
    lagged.logxs <- rep(1,p)
    lagged.psis <- rep(1,q)
    
    x <-rep(NA,n)	# we save last n entries of xs into x 
    psi <-rep(NA,n)	# we save last n entries of psis into psi
    
    
    ########################################################################
    ## Generate errors eps from a positive-valued distribution
    ## Choices are: exponential, Weibull, or Gamma
    ########################################################################
    
    if (feps=="exponential") 
    {
      ## Generate errors from exponential distribution
      eps <- eps.exp <- rexp(nt,par1)
    } else if (feps=="weibull") {		
      ## Generate errors from weibull distribution			
      eps <- eps.weib <- rweibull(nt, shape=par1, scale=par2) 	 
    } else if (feps=="gamma") {
      ## Generate errors from gamma distribution			
      # randomly generate errors from gamma distribution			
      eps <- eps.gamma <- rgamma(nt,shape=par1,scale=par2)
    }
    
    
    ########################################################################
    ## Compute psis and xs
    ########################################################################
    
    for (t in 1:mxpq){
      psis[t] = omega
      logxs[t]= psis[t] + log(eps[t]) -log(mean(eps))	
      xs[t]=exp(logxs[t])
    }
    
    for (t in (mxpq+1):nt) 
    {
      lagged.logxs <- c(logxs[(t-1):(t-p)])
      lagged.psis <- c(psis[(t-1):(t-q)])
      psis[t] = omega + alpha%*%lagged.logxs + beta%*%lagged.psis
      logxs[t] = psis[t] + log(eps[t]) 	-log(mean(eps))			
      xs[t] = exp(logxs[t] )
    }
    
  }
  ## End Durations Simulation of psis and xs
  
  ## Save last n values as time series into x and psi (discard first nb as burn-in)		
  x <- ts(xs[(nb+1):nt])	 	
  psi <- ts(psis[(nb+1):nt])	
  
  ## output as dataframe, Col 1 has x and Col 2 has psi		
  output <- as.data.frame(cbind(x,psi))
}
