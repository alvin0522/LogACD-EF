setwd("path where all codes are stored")
source("core_funs.R")
source("user_define.R")
source("AVRE.R")
source("AISRE.R")
source("NESE.R")
source("LogAcd_EF.R")

if("R.matlab" %in% rownames(installed.packages()) == FALSE)
  install.packages("R.matlab")
library(R.matlab)
Matlab$startServer( )        ### Open Matlab Server                                 
matlab <- Matlab()          
isOpen <- open(matlab)       ### Connect R and Matlab
isOpen  


########### Non-Penalized EF Approach ############
## set up function parameters ##

### Five scenarios in simulation study
#omega = 0.25, alpha = 0.05    Gamma(0.6, 0.7)
#omega = 0.04, alpha = 0.05, beta = 0.75)  Exp(1)
#omega = 1, alpha = 0.05, beta = 0.75)   Weibull(0.4, 0.5)
#omega = 0.5, alpha1 = 0.05, alpha2 = 0.1, beta = 0.6)  Weibull(0.9, 0.9)
#omega = 0.15, alpha1 = 0.05, alpha2 = 0.1, beta1 = -0.05, beta2 = 0.7)  Gamma(0.5, 0.8)

omega <- 0.25
alpha <- 0.05
beta <- 0
sim_dist <- "weibull"  # 'gamma', 'weibull'. 'define par1=1 and sim_dist='gamma' to get exp(par2) #
par1 <-  0.9
par2 <-  0.9
nsim <-  2
len_sim <- 7500
nburn <- 500
method <- "NESE"  ## Three methods: "NESE' 'AVRE' 'AISRE'
n_interation <- 50

set.seed(123457)

## call functions - Reproduce Raw Estimates for Figure 8 ##
myres <- LogAcd_EF(omega = omega, alpha = alpha, beta = beta, 
                  sim_dist = sim_dist, par1 = par1, par2 = par2,
                  nsim = nsim, len_sim = len_sim, nburn = nburn,
                  method = method, n_interation = n_interation)

res_summary <- t(apply(X = myres$estimation, MARGIN = 2, FUN = quantile, 
                       probs = c(0.05,0.25, 0.5, 0.75,0.95), na.rm=T))


########### Penalized EF Approach ############
omega <- 0.25
alpha <- c(0.2, 0.10)
beta <- 0
sim_dist <- "gamma"  # 'gamma', 'weibull'. 'define par1=1 and sim_dist='gamma' to get exp(par2) #
par1 <-  0.5
par2 <-  0.6
nsim <-  2
len_sim <- 7500
nburn <- 500
method <- "AVRE"  ## Two methods available: "NESE' 'AVRE' 
n_interation <- 50
lambda <- 0.3
a <- 3.7
p <- 20
q <- 0
cutoff <- 10^-3

set.seed(123457)
myres <- LogAcd_EF_Pen(omega = omega, alpha = alpha, beta = beta, lambda= lambda, a=a,
                          p = p, q= q, sim_dist = sim_dist, par1 = par1, par2 = par2,
                          nsim = nsim,len_sim = len_sim, nburn = nburn, method = method, 
                          n_interation = n_interation, cutoff = cutoff)


