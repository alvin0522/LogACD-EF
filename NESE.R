## Method Nonlinear Equation Solver Estimation ##
# Inputs are:
# x: durations (vector)
# initval: initial values for Log ACD model
# p: LogACD model, first parameter
# q: LogACD model, second parameter
# momente: first four moments (vector)

# Outputs are:
# estimation (vector)


nleqn.matlab <- function(x, initval, p, q, momente){
  ini.para=c(as.numeric(initval),0.01)                              
  x.write = as.numeric(c(x,p,q))
  pdim <- p+q+2
  dist.write <- momente
  write(x.write, "data_info.csv", ncolumns=1,sep=",")
  write(dist.write, "moments.csv" ,ncolumns=1, sep=",")
  
  ## Feed R variables to Matlab
  #wd <- getwd()
  #setVariable(matlab, wd=wd)
  setVariable(matlab, para=ini.para, p=p, q=q)
  
  
  evaluate(matlab,  
           "npara=p+q+2;",
           "lb=repmat(-1/0,1,npara);",
           "lb(2:npara-1)=-1;",
           "ub=repmat(1/0,1,npara);",
           "ub(2:npara)=1;",
           "lb(npara)=0;",
           "y=lsqnonlin(@nese,para,lb,ub);"
  )
  
  res <- getVariable(matlab,"y")                                       
  result <- as.numeric(res$y)[-pdim]
  return(result)
}



## Method Nonlinear Equation Solver Estimation - Penalized Version##
# Inputs are:
# x: durations (vector)
# initval: initial values for Log ACD model
# p: LogACD model, first parameter
# q: LogACD model, second parameter
# momente: first four moments (vector)

# Outputs are:
# estimation (vector)


nleqn.matlab.pen <- function(x, initval, p, q, momente, l1, a){
  ini.para=c(as.numeric(initval),0.01)                              
  x.write = as.numeric(c(x,p,q))
  pdim <- p+q+2
  dist.write <- momente
  lambda <- c(l1, a)
  write(x.write, "data_info1.csv", ncolumns=1,sep=",")
  write(dist.write, "moments1.csv" ,ncolumns=1, sep=",")
  write(lambda, "lambda1.csv", ncolumns=1,sep=",")
  write(as.numeric(initval), "initial1.csv", ncolumns=1,sep=",")
  ## Feed R variables to Matlab
  wd <- getwd()
  setVariable(matlab, wd=wd)
  setVariable(matlab, para=ini.para, p=p, q=q)
  
  
  evaluate(matlab,  
           "npara=p+q+2;",
           "lb=repmat(-1/0,1,npara);",
           "lb(2:npara-1)=-1;",
           "ub=repmat(1/0,1,npara);",
           "ub(2:npara)=1;",
           "lb(npara)=0;",
           "y=lsqnonlin(@nese_pen,para,lb,ub);"
  )
  
  res <- getVariable(matlab,"y")                                       
  result <- as.numeric(res$y)[-pdim]
  return(result)
}
