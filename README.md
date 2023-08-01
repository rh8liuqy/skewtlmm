# Skew-t Linear Mixed Model fit with Distributed AECM Algorithm
In order to use the script, you will need to download the STDEC_script.R and fastLik.cpp files, along with updating your packages with `install.packages(c('parallel', 'RcppArmadillo', 'Rcpp', 'purrr', 'dplyr', 'mvtnorm'))` . 
# Use
The primary script, utilized to fit the skew- $t$ linear mixed model is uploaded here. This function takes primary arguments of $Y$, a list of response vectors per subject with overall list length $N$, $X$, a list of covariate matrices per subject with list length $N$, $T_i$, a list of observation times per subject, and $n_i$, a list of number of observations per subject. Initial values should be provided to the function of all parameters in the model. As the $AECM$ algorithm only guarantees convergence to a local minima of the objective function, the script should ideally be rerun multiple times with different starting values in order to ensure convergence to the global optima. In practice, users may want to utilize smaller tolerances to ensure convergence. 

# Install the package

```r
#Install the package.
devtools::install_github("rh8liuqy/skewtlmm")
```

# Example


```r
library(skewtlmm)
#Set a seed
set.seed(1)
## Load library
library(purrr)
library(mvtnorm)
library(dplyr)
library(parallel)
library(doParallel)
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
    bi = split(rst(N,xi=0,s=1,la=3,df=50),1:N)
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

#The below code utilizes 12 cores.
output <- distributed_EM_model(numCores = 12,
                               nwaitf = 3,
                               Yf = Y,
                               Xf = X,
                               Zf = Z,
                               nif = ni,
                               T_if = T_i,
                               option = "skewT-DEC")

output <- distributed_EM_model(numCores = 12,
                               nwaitf = 3,
                               Yf = Y,
                               Xf = X,
                               Zf = Z,
                               nif = ni,
                               T_if = T_i,
                               option = "skewN-DEC")

output <- distributed_EM_model(numCores = 12,
                               nwaitf = 3,
                               Yf = Y,
                               Xf = X,
                               Zf = Z,
                               nif = ni,
                               T_if = T_i,
                               option = "normal-DEC")

output <- distributed_EM_model(numCores = 12,
                               nwaitf = 3,
                               Yf = Y,
                               Xf = X,
                               Zf = Z,
                               nif = ni,
                               T_if = T_i,
                               option = "skewT-normal")

output <- distributed_EM_model(numCores = 12,
                               nwaitf = 3,
                               Yf = Y,
                               Xf = X,
                               Zf = Z,
                               nif = ni,
                               T_if = T_i,
                               option = "skewN-normal")

output <- distributed_EM_model(numCores = 12,
                               nwaitf = 3,
                               Yf = Y,
                               Xf = X,
                               Zf = Z,
                               nif = ni,
                               T_if = T_i,
                               option = "normal-normal")

print("all success!")
```
