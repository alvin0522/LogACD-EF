setwd("path where all codes are stored")
source("core_funs.R")
source("user_define.R")
source("recursive_mat.R")
source("recursive_scalar.R")
source("nleqn.R")
source("LogAcd_EF.R")

if("R.matlab" %in% rownames(installed.packages()) == FALSE)
  install.packages("R.matlab")
library(R.matlab)
Matlab$startServer( )        ### Open Matlab Server                                 
matlab <- Matlab()          
isOpen <- open(matlab)       ### Connect R and Matlab
isOpen  

## set up function parameters ##
omega <- 0.15
alpha <- c(0.1, -0.05)
beta  <- c(0.05, 0.7)
sim_dist <- "gamma"  # 'gamma', 'weibull'. 'define par1=1 and sim_dist='gamma' to get exp(par2) #
par1 <-  0.5
par2 <-  0.8
nsim <-  250
len_sim <- 7500
nburn <- 500
method <- "AVRE"  ## Three methods: "NESE' 'AVRE' 'AISRE'
n_interation <- 50

set.seed(123457)

## call functions - Reproduce Table 1 ##
myres <- LogAcd_EF(omega = omega, alpha = alpha, beta = beta, 
                  sim_dist = sim_dist, par1 = par1, par2 = par2,
                  nsim = nsim, len_sim = len_sim, nburn = nburn,
                  method = method, n_interation = n_interation)

res_summary <- t(apply(X = myres$estimation, MARGIN = 2, FUN = quantile, 
                       probs = c(0.05,0.25, 0.5, 0.75,0.95), na.rm=T))
