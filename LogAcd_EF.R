## Log ACD(p,q) model -- Estimating Functions Approach
## Inputs are:
## omega, alpha, beta: parameters described in Log ACD(p,q) model
## sim_dist: distribution of epsilon
## par1, par2: parameters of epsilon distribution
## nsim: number of simulations
## len_sim: length of simulated durations
## nburn: number of burn-in durations
## method: method to estimate parameters. Three methods- "NESE", "AVRE", "AISRE"
## n_interation: maximum number of iteration, if "AISRE" is specified

## Outputs are:
## parameter estiamtion (vector)


LogAcd_EF <- function(omega, alpha, beta, sim_dist, par1 = NULL,
                      par2 = NULL, nsim, len_sim, nburn, method, n_interation=NULL){
  p <- length(alpha)
  q <- length(beta)
  if(sum(as.numeric(beta == 0)) == 1){
    q <-0
  }else {
    q <- length(beta)
  }
  pdim <- p+q+1
  EF_para <- matrix(0, nrow=nsim, ncol=pdim)
  EF_se <- matrix(0, nrow=nsim, ncol=pdim)
  
  if(method == "NESE"){
    wd <- getwd()
    setVariable(matlab, wd=wd)
    evaluate(matlab, "cd(wd);")
  }
  
  for(i in 1:nsim){
    dur <- fsim.logacd(len_sim, nburn, omega, alpha, beta,
                        feps = sim_dist, par1 = par1, par2=par2)
    x <- dur$x
    est <- finitval(x,p,q)
    if (q==0) est <- est[-(pdim+1)]
    psih <- psi_cal(x,est,p,q)
    eps <- x/exp(psih)
    ed <- eps_dist2(eps)
    momente <- c(ed$mue, ed$vare, ed$skewe, ed$kurte)
    difference <- log(momente[1]) - mean(log(eps))
    
    if(method=="AVRE"){
      
      EF_est <- recursive_mat(x,est,p,q,momente)
      if(q==0) {
        EF_est[1] <- EF_est[1] + difference
      }
      if(q >0) {
        EF_est[1] <- EF_est[1] + (1-sum(EF_est[(p+2):pdim]))*difference
      }
      
    }else if(method == "AISRE"){
      
      EF_est <- festeq.logacd1.online.pq(x, est,p,q,momente,iteration=n_interation)  
      if(q==0) {
        EF_est[1] <- EF_est[1] + difference
      }
      if(q >0) {
        EF_est[1] <- EF_est[1] + (1-sum(EF_est[(p+2):pdim]))*difference
      }
      
    }else {
      EF_est <- nleqn.matlab(x, est, p, q, momente)     
    }
    
    se <- StdErrCal(x = x, initval = EF_est, p = p, q = q, momente)
    EF_para[i,] <- EF_est
    EF_se[i,] <- se
    #print(i)
    
  }
  
  return(list(estimation = EF_para, se=EF_se))
  
}


########### Penalized EF Approach #############
LogAcd_EF_Pen <- function(omega, alpha, beta, lambda, a, p, q, sim_dist, par1 = NULL,
                          par2 = NULL, nsim, len_sim, nburn, method, 
                          n_interation=NULL, cutoff){
  if(sum(as.numeric(beta == 0)) == 1){
    q <-0
  }else {
    q <- length(beta)
  }
  pdim <- p+q+1
  EF_para <- matrix(0, nrow=nsim, ncol=pdim)
  EF_se <- matrix(0, nrow=nsim, ncol=pdim)
  
  if(method == "NESE"){
    wd <- getwd()
    setVariable(matlab, wd=wd)
    evaluate(matlab, "cd(wd);")
  }
  
  for(i in 1:nsim){
    dur <- fsim.logacd(len_sim, nburn, omega, alpha, beta,
                       feps = sim_dist, par1 = par1, par2=par2)
    x <- dur$x
    est <- finitval(x,p,q)
    if (q==0) est <- est[-(pdim+1)]
    psih <- psi_cal(x,est,p,q)
    eps <- x/exp(psih)
    ed <- eps_dist2(eps)
    momente <- c(ed$mue, ed$vare, ed$skewe, ed$kurte)
    difference <- log(momente[1]) - mean(log(eps))
    
    if(method=="AVRE"){
      
      EF_est <- recursive_mat2(x, est, p, q, momente, lambda, a)
      if(q==0) {
        EF_est[1] <- EF_est[1] + difference
      }
      if(q >0) {
        EF_est[1] <- EF_est[1] + (1-sum(EF_est[(p+2):pdim]))*difference
      }
      
    }else {
      EF_est <- nleqn.matlab.pen(x, est, p, q, momente, lambda, a)
    }
    EF_est <- ifelse(abs(EF_est) < cutoff, 0, EF_est)
    se <- StdErrCal(x = x, initval = EF_est, p = p, q = q, momente)
    EF_para[i,] <- EF_est
    EF_se[i,] <- se
    #print(i)
    
  }
  
  
  return(list(estimation = EF_para, se=EF_se))
  
}

