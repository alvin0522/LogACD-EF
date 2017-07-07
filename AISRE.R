## Method Approximate Interated Scalar Recursive Estimation (AISRE) ##
# Inputs are:
# x: durations (vector)
# initval: initial values for Log ACD model
# p: LogACD model, first parameter
# q: LogACD model, second parameter
# momente: first four moments (vector)
# iteration: Maximum value 

# Outputs are:
# estimation (vector)

festeq.logacd1.online.pq <- function(x,initval,p,q,momente,iteration){
  
  n <- length(x)
  pdim <- 1+p+q   # number of parameters to be estimated
  mxpq <- max(p,q)
  
  mue <- momente[1]
  vare <- momente[2]
  skewe <- momente[3]
  kurte <- momente[4]
  
  ########################################################################    
  ## INITIALIZATION 
  ########################################################################
  
  ## psi hat
  psih <- rep(1,n)
  
  ## Parameter estimates for each iteration
  initial<- as.numeric(initval)
  thehat <- matrix(rep(initial,n),ncol=n,nrow=pdim,byrow = F)
  
  
  
  ########################################################################
  ## RECURSIVE ONLINE ESTIMATION ON THETA
  ########################################################################
  ## Put initial values into initial positions of arrays
  thehat[,1]=initial
  kmat <- diag(solve(InitMat(x,p,q)))
  kmatold <- rep(1,pdim)
  theta <- matrix(NA, ncol=iteration, nrow=pdim)
  err <- rep(1,pdim)
  
  for(iter in 1:iteration){
    
    #psih[1]=thehat[1,1]    # omega   
    for( j in 1:pdim  ){
      
      for (t in 1:mxpq){
        psih[t] = thehat[1,1]
      }
      
      
      for( t in (mxpq+1):n ){
        omega <- thehat[1,t-1]
        alpha <- as.matrix(thehat[2:(p+1),(t-1)],nrow=p)
        if(q !=0){
          beta <- as.matrix(thehat[(p+2):pdim,t-1],nrow=q)
        }
          pos1 <- t-1
          pos2 <- t-p
          pos3 <- t-q
          
        psih[t]=psit(x,psih,omega, alpha, beta, pos1, pos2, pos3, q)
			
        ## Define psi								
        
        kmatold[j] <- kmat[j]
        ## derpsi w.r.t alpha ##
        derpsi.v <- c(1,log(x[(t-1):(t-p)]),psih[(t-1):(t-q)])
        derpsi <- derpsi.v[j]
        der2psi=0
        
        
        ## mu(t), sigsq(t), gamma(t), kappa(t)							
        mu=mue*exp(psih[t]) 							
        sigsq=vare*exp(2*psih[t]) 							
        gamma=skewe*exp(3*psih[t])  							
        kappa=kurte*exp(4*psih[t])  
        
        
        ## First Derivatives of mu(t),sigsq(t), gamma(t) and kappa(t)												
        dermu=mue*exp(psih[t])*derpsi							
        dersigsq=2*vare*exp(2*psih[t])*derpsi
        dergamma=3*skewe*exp(3*psih[t])*derpsi
        derkappa=4*kurte*exp(4*psih[t])*derpsi
        
        ## Second Derivatives of mu(t) and sigsq(t)  						
        #der2mu=mue*der2psi							
        #der2sigsq=2*vare*(der2psi*psih[t]+derpsi**2)
        
        der2mu <- mu*((derpsi)%*%t(derpsi) + der2psi)
        der2sigsq <- 2*sigsq*(der2psi + 2*(derpsi)%*%t(derpsi))
        
        ## END CHANGE_mod
        ##########################################################
        
        
        ## Compute m(t) and M(t)							
        m = x[t]-mu       							
        qm = m**2-sigsq
        
        ## Compute Quadratic variations of m(t) and M(t) and 
        ## covariance of (m(t), M(t))							
        vm = sigsq
        #vqm = (vm**2)*(kappa-1)  #Note: kappa is NOT excess kurtosis
        vqm=kappa-vm**2
        #vmqm = (sigsq**(3/2))*gamma
        vmqm=gamma
        
        
        ## Define rho^2(t) 															
        termr = 1-((vmqm**2)/(vm*vqm))								
        rho = 1/termr
        
        ## Define eta 								
        eta = vmqm/(vm*vqm)		
        
        ## Define vectors astr and bstr						
        astr = rho*(-dermu/vm + dersigsq*eta)							
        bstr = rho*(dermu*eta - dersigsq/vqm)
        
        ## derivatives of m and qm and their variations
        derm = -dermu
        derqm = 2*m*derm-dersigsq
        
        dervm = dersigsq
        #dervqm = 2*(kappa-1)*sigsq*dersigsq+(sigsq*2)*derkappa
        dervqm=derkappa - 2*sigsq*dersigsq
        #dervmqm = (sigsq**(3/2))*dergamma+(3/2)*gamma*sqrt(sigsq)*dersigsq
        dervmqm = dergamma
        
        ## Derivatives of eta and rho
        anum = vm*vqm*dervmqm-vmqm*(vm*dervqm+vqm*dervm)
        denom = (vm**2)*(vqm**2)
        dereta = anum/denom
        
        v = vm*vqm -vmqm**2
        u = vm*vqm
        dv = vm*dervqm+vqm*dervm-2*vmqm*dervmqm
        du = vm*dervqm+vqm*dervm
        derrho = (v*du-u*dv)/(v**2)
        
        ## Derivatives of astr and bstr  						
        ## For astr						
        terma1 = (vm*der2mu-dermu*dervm)/(vm**2)
        terma2 = dersigsq*dereta+eta*der2sigsq
        terma3 = (-dermu/vm + dersigsq*eta)*derrho						
        derastr = -rho*terma1 + rho*terma2 + terma3
        
        
        ## For bstr						
        termb1 = dermu*dereta+eta*der2mu												
        termb2 = (vqm*der2sigsq-dersigsq*dervqm)/(vqm**2)
        termb3 = (dermu*eta-dersigsq/vqm)*derrho
        derbstr = rho*termb1 - rho*termb2 + termb3
        
        ##
        aterm = derastr*m+derm*astr+derbstr*qm+derqm*bstr
        kmat[j] = kmatold[j]/(1-aterm*kmatold[j])
        
        ##
        ## Add numerical control
        if(is.na(kmat[j])) kmat[j] <- kmatold[j]
        
        ## compute new theta						
        termt1 = astr*m + bstr*qm 
        temp <- thehat[j,t-1] + kmat[j]*termt1
        if((1-sum(temp[-1])) < 0.001) thehat[j,t] <- thehat[j,t-1]
        else thehat[j,t]  = temp
        
        
      }
      
    }
    theta[,iter] <- thehat[,n]
    for(m in 1:pdim){
      if(is.na(theta[m,iter])) {
        flag <- 0
        for ( w in 1: iter-1){
          if(flag == 0 & !is.na(theta[m,iter-w])){
            theta[m,iter] <- theta[m,iter-w]
            flag <- 1
          }
          else theta[m,iter] <- initial[m]
        }
      }
    }
    if(iter >=2) err <- abs(theta[,iter]-theta[,iter-1])
    change <- sum(err < rep(0.001, pdim))
    if(change == pdim) break
  }
  index <- apply(thehat, MARGIN = 2, is.na)
  ind <- apply(index, MARGIN = 2, sum)
  ind_f <- ind == 0
  theta_est <- matrix(thehat[,ind_f],nrow=pdim)
  nc <- ncol(theta_est)
  #return(est=theta_est[,nc])
  return(est=theta_est[,nc])
  
}

#### LogACD2 Model ###
festeq.logacd2.online.pq <- function(x,initval,p,q,mue,vare,skewe,kurte,iteration){
  
  n <- length(x)
  pdim <- 1+p+q   # number of parameters to be estimated
  mxpq <- max(p,q)
  
  ########################################################################    
  ## INITIALIZATION 
  ########################################################################
  
  ## psi hat
  psih <- rep(1,n)
  
  ## Parameter estimates for each iteration
  initial<- as.numeric(initval)
  thehat <- matrix(rep(initial,n),ncol=n,nrow=pdim,byrow = F)
  
  
  
  ########################################################################
  ## RECURSIVE ONLINE ESTIMATION ON THETA
  ########################################################################
  ## Put initial values into initial positions of arrays
  thehat[,1]=initial
  kmat <- diag(solve(InitMat2(x,p,q,difference)))
  kmatold <- rep(1,pdim)
  theta <- matrix(NA, ncol=iteration, nrow=pdim)
  err <- rep(1,pdim)
  
  for(iter in 1:iteration){
    
    #psih[1]=thehat[1,1]    # omega   
    for( j in 1:pdim  ){
      
      for (t in 1:mxpq){
        psih[t] = thehat[1,1]
      }
      
      
      for( t in (mxpq+1):n ){
        lagged.psih <- as.matrix(c(psih[(t-1):(t-q)]))
        eps <- x/exp(psih)
        alpha.term <- as.matrix(eps[(t-1):(t-p)])
        omega <- thehat[1,t-1]
        alpha <- as.matrix(thehat[2:(p+1),(t-1)],nrow=p)
        beta <- as.matrix(thehat[(p+2):pdim,t-1],nrow=q)
        psih[t]=omega+t(alpha)%*%alpha.term+t(beta)%*%lagged.psih
        ## CHANGE_mod (change psi, derivative of psi and second derivative of psi)
        ## change until END CHANGE_mod      				
        ## Define psi								
        
        kmatold[j] <- kmat[j]
        ## derpsi w.r.t alpha ##
        
        derpsi.v <- c(1,eps[(t-1):(t-p)],psih[(t-1):(t-q)])
        derpsi <- derpsi.v[j]
        der2psi=0
        
        
        ## mu(t), sigsq(t), gamma(t), kappa(t)							
        mu=exp(psih[t]) 							
        sigsq=vare*exp(2*psih[t])/mue^2 							
        gamma=skewe*exp(3*psih[t])/mue^3  # recall skewe is third central moment							
        kappa=kurte*exp(4*psih[t])/mue^4  # recall kurte is fourth central moment
        
        
        ## First Derivatives of mu(t),sigsq(t), gamma(t) and kappa(t)												
        dermu=exp(psih[t])*derpsi							
        dersigsq=2*vare*exp(2*psih[t])*derpsi/mue^2
        dergamma=3*skewe*exp(3*psih[t])*derpsi/mue^3
        derkappa=4*kurte*exp(4*psih[t])*derpsi/mue^4
        
        ## Second Derivatives of mu(t) and sigsq(t)  						
        #der2mu=mue*der2psi							
        #der2sigsq=2*vare*(der2psi*psih[t]+derpsi**2)
        
        der2mu <- mu*((derpsi)%*%t(derpsi) + der2psi)
        der2sigsq <- 2*sigsq*(der2psi + 2*(derpsi)%*%t(derpsi))
        
        ## END CHANGE_mod
        ##########################################################
        
        
        ## Compute m(t) and M(t)							
        m = x[t]-mu       							
        qm = m**2-sigsq
        
        ## Compute Quadratic variations of m(t) and M(t) and 
        ## covariance of (m(t), M(t))							
        vm = sigsq
        #vqm = (vm**2)*(kappa-1)  #Note: kappa is NOT excess kurtosis
        vqm=kappa-vm**2
        #vmqm = (sigsq**(3/2))*gamma
        vmqm=gamma
        
        
        ## Define rho^2(t) 															
        termr = 1-((vmqm**2)/(vm*vqm))								
        rho = 1/termr
        
        ## Define eta 								
        eta = vmqm/(vm*vqm)		
        
        ## Define vectors astr and bstr						
        astr = rho*(-dermu/vm + dersigsq*eta)							
        bstr = rho*(dermu*eta - dersigsq/vqm)
        
        ## derivatives of m and qm and their variations
        derm = -dermu
        derqm = 2*m*derm-dersigsq
        
        dervm = dersigsq
        #dervqm = 2*(kappa-1)*sigsq*dersigsq+(sigsq*2)*derkappa
        dervqm=derkappa - 2*sigsq*dersigsq
        #dervmqm = (sigsq**(3/2))*dergamma+(3/2)*gamma*sqrt(sigsq)*dersigsq
        dervmqm = dergamma
        
        ## Derivatives of eta and rho
        anum = vm*vqm*dervmqm-vmqm*(vm*dervqm+vqm*dervm)
        denom = (vm**2)*(vqm**2)
        dereta = anum/denom
        
        v = vm*vqm -vmqm**2
        u = vm*vqm
        dv = vm*dervqm+vqm*dervm-2*vmqm*dervmqm
        du = vm*dervqm+vqm*dervm
        derrho = (v*du-u*dv)/(v**2)
        
        ## Derivatives of astr and bstr  						
        ## For astr						
        terma1 = (vm*der2mu-dermu*dervm)/(vm**2)
        terma2 = dersigsq*dereta+eta*der2sigsq
        terma3 = (-dermu/vm + dersigsq*eta)*derrho						
        derastr = -rho*terma1 + rho*terma2 + terma3
        
        
        ## For bstr						
        termb1 = dermu*dereta+eta*der2mu												
        termb2 = (vqm*der2sigsq-dersigsq*dervqm)/(vqm**2)
        termb3 = (dermu*eta-dersigsq/vqm)*derrho
        derbstr = rho*termb1 - rho*termb2 + termb3
        
        ##
        aterm = derastr*m+derm*astr+derbstr*qm+derqm*bstr
        kmat[j] = kmatold[j]/(1-aterm*kmatold[j])
        
        ##
        ## Add numerical control
        if(is.na(kmat[j])) kmat[j] <- kmatold[j]
        
        ## compute new theta						
        termt1 = astr*m + bstr*qm 
        temp <- thehat[j,t-1] + kmat[j]*termt1
        if((1-sum(temp[-1])) < 0.001) thehat[j,t] <- thehat[j,t-1]
        else thehat[j,t]  = temp
        
        
      }
      
    }
    theta[,iter] <- thehat[,n]
    for(m in 1:pdim){
      if(is.na(theta[m,iter])) {
        flag <- 0
        for ( w in 1: iter-1){
          if(flag == 0 & !is.na(theta[m,iter-w])){
            theta[m,iter] <- theta[m,iter-w]
            flag <- 1
          }
          else theta[m,iter] <- initial[m]
        }
      }
    }
    if(iter >=2) err <- abs(theta[,iter]-theta[,iter-1])
    change <- sum(err < rep(0.001, pdim))
    if(change == pdim) break
  }
  index <- apply(thehat, MARGIN = 2, is.na)
  ind <- apply(index, MARGIN = 2, sum)
  ind_f <- ind == 0
  theta_est <- matrix(thehat[,ind_f],nrow=pdim)
  nc <- ncol(theta_est)
  #return(est=theta_est[,nc])
  return(est=theta_est[,nc])
  
}
