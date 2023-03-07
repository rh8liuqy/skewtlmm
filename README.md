# Skew-t Linear Mixed Model fit with Distributed AECM Algorithm
In order to use the script, you will need to download the STDEC_script.R and fastLik.cpp files, along with updating your packages with `install.packages(c('parallel', 'RcppArmadillo', 'Rcpp', 'purrr', 'dplyr', 'mvtnorm'))` . 
# Use
The primary script, utilized to fit the skew- $t$ linear mixed model is uploaded here. This function takes primary arguments of $Y$, a list of response vectors per subject with overall list length $N$, $X$, a list of covariate matrices per subject with list length $N$, $T_i$, a list of observation times per subject, and $n_i$, a list of number of observations per subject. Initial values should be provided to the function of all parameters in the model. As the $AECM$ algorithm only guarantees convergence to a local minima of the objective function, the script should ideally be rerun multiple times with different starting values in order to ensure convergence to the global optima. In practice, users may want to utilize smaller tolerances to ensure convergence. 

# Example


```r
#Set a seed
set.seed(1352)
## Load library
library(purrr)
library(mvtnorm)
library(dplyr)
#Function to generate skew-normal deviates
rsn = function(n,xi=0,s=1,la)
{
  del = la/sqrt(1+la^2)
  z1 = abs(rnorm(n))
  z2 = rnorm(n)
  xi+s*(del*z1+sqrt(1-del^2)*z2)
}

#Function to generate skew-t deviates
rst = function(n,xi=0,s=1,la,df)
{
  xi + rsn(n,la=la)*s/sqrt(rgamma(n,df/2,df/2))
}

#Data-generating function 
gen_DS = function(N = 1e5) {
  bi = split(rst(N,xi=0,s=1,la=3,df=4),1:N)
  ni = rpois(N, 5)+2
  tList = lapply(ni, function(x) cumsum(rpois(x,3)+1))
  covar_params = c(0.5,0.9)
  ei = pmap(list(ni,tList), function(n,t_i) rmvt(n = 1, sigma = covar_params[1]^abs(outer(t_i,t_i,'-'))^covar_params[2],df = 4) %>% t)
  
  xList = lapply(tList, function(x) matrix(c(rep(1, length(x)),rnorm(length(x)), rnorm(length(x))), ncol = 3))
  yList = pmap(list(ei, bi, xList), function(e,b, xl) -2 + 2*xl[,2] - 2*xl[,3]+ as.numeric(b) + e)
  zList = lapply(yList, function(x) matrix(c(rep(1, length(x))),ncol = 1))
  return(list(yList, xList, zList,tList,ni))
}

#Create simulated data
simulated_data = gen_DS(N = 1e5)
#Setup data values
Y = simulated_data[[1]]
X = simulated_data[[2]]
Z = simulated_data[[3]]
T_i = simulated_data[[4]]
ni = simulated_data[[5]]
N = length(Y)

#The below code utilizes 12 cores to fit the model with 3 different fractional return hyperparameters: 0.25, 0.5, and 0.75. 

library(purrr)
library(parallel)
library(mvtnorm)
library(doParallel)

#setup parallel cluster
cl <- makeCluster(numCores)
setDefaultCluster(cl)
registerDoParallel(cl)
idList = split(1:N, cut(1:N, c(0,sample(1:N, 12-1),N)))
clusterExport(cl, c('idList', 'pmap'))

yL = Y
xL = X
zL = Z 
niL = ni 
tiL = T_i 

#export data to the cluster we just set up
clusterExport(cl, c('yL', 'xL', 'zL', 'niL', 'tiL'))

clusterApply(cl, 1:numCores, function(x){
  ind_L = idList[[x]]
  Y <<- yL[ind_L]
  X <<- xL[ind_L]
  Z <<- zL[ind_L]
  ni <<- niL[ind_L]
  T_i <<- tiL[ind_L]
  return(TRUE)
})

numCores = 12
nwait = 3
source('STDEC_script.R')
gamma_25 = para

nwait = 6
source('STDEC_script.R')
gamma_50 = para

nwait = 9
source('STDEC_script.R')
gamma_75 = para

#Time
print(c(gamma_25$time[1], gamma_50$time[1], gamma_75$time[1]))
#Iterations
print(c(gamma_25$time[2], gamma_50$time[2], gamma_75$time[2]))
            

```
