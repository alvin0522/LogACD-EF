# LogACD-EF
==========
##Overview
The code allows the user to carry out Log ACD$(p,q)$ modeling of time series of durations using martinglae estimating functions, as discussed in our manuscript. Specifically, this code will enable the user to generate time series of durations as discussed in section 7 and to reproduce the results in  Table 1. With slight modification, the code is used for the real data analysis shown in section 8. 

##Required Software
*  R/Rstudio 
*  Matlab 

##Required R Packages
Please install the following packages before you execute the sample code

* R.matlab
* moments
* MASS

##How to Use

1. Place all seven R files and one matlab file in the same folder.
2. Add path-to-folder to Matlab.
3. Open example_code.R and change your working directory:
```{evaluate = FALSE}
setwd("path where all codes are stored")
```
4. Source all necessary codes.
5. Connect R and Matlab:
```{r load_packages}
library(R.matlab)
Matlab$startServer( )        ### Open Matlab Server                                 
matlab <- Matlab()          
isOpen <- open(matlab)       ### Connect R and Matlab
isOpen                       ### If successfully connected, TRUE will be returned
```
6. Set up simulation parameters:

  Arguments  | Description
  -----------|------------
  omega, alpha, beta | parameters of the Log ACD$(p,q)$ model
  sim_dist | distribution of errors epsilon
  par1, par2 | parameters of the epsilon distribution
  nsim | number of simulations
  len_sim | length of simulated durations
  nburn | number of burn-in durations
  method | method to estimate parameters ("NESE", "AVRE", "AISRE")
  n_interation | maximum number of iterations, if "AISRE" is specified

7. To reproduce our results, set seed as $123457$

