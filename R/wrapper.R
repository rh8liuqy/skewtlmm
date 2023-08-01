library(parallel)
library(RcppArmadillo)
library(Rcpp)

distributed_EM_model <- function(numCores,
                                 nwaitf,
                                 Yf,
                                 Xf,
                                 Zf,
                                 nif,
                                 T_if,
                                 option) {
    if (numCores < nwaitf) {
        stop("Error: numCores must be larger than or equal to nwaitf!")
        return(NULL)
    }
    N <<- length(Y)
    numCores <<- numCores

    #setup parallel cluster
    cl <<- makeCluster(numCores)
    setDefaultCluster(cl)
    registerDoParallel(cl)
    idList <<- split(1:N, cut(1:N, c(0,sample(1:N, numCores-1),N)))
    clusterExport(cl, c('idList', 'pmap'))
    clusterEvalQ(cl, {library(purrr)
        library(parallel)
        library(mvtnorm)
        library(doParallel)})
    Y <<- Yf
    X <<- Xf
    Z <<- Zf
    ni <<- nif
    T_i <<- T_if
    yL <<- Y
    xL <<- X
    zL <<- Z
    niL <<- ni
    tiL <<- T_i

    p <- ncol(X[[1]])
    be_init <- numeric(p)
    be_init[1] <- 1.0

    #export data to the cluster we just set up
    clusterExport(cl, c('yL', 'xL', 'zL', 'niL', 'tiL'))

    clusterApply(cl, 1:numCores, function(x){
        ind_L <<- idList[[x]]
        Y <<- yL[ind_L]
        X <<- xL[ind_L]
        Z <<- zL[ind_L]
        ni <<- niL[ind_L]
        T_i <<- tiL[ind_L]
        return(TRUE)
    })

    nwait <<- nwaitf
    if (option == "skewT-DEC") {
        script_path <- system.file("/Reuben_Retnam_Code/STDEC_script.R",package = "skewtlmm")
        source(script_path)
    }
    else if (option == "skewT-normal") {
        script_path <- system.file("/Reuben_Retnam_Code/STNormal_script.R",package = "skewtlmm")
        source(script_path)
        para <- STLMM.AECM.ND(Y=Y,
                              X=X,
                              Z=Z,
                              N=N,
                              ni=ni,
                              T_i=T_i,
                              be=be_init,
                              la=c(1),
                              s2=2,
                              Ga=matrix(c(0.35),1,1),
                              nu=7,
                              tol=1e-4,
                              max.iter=5e5,
                              per=1)
    }
    else if (option == "skewN-DEC") {
        script_path <- system.file("/Reuben_Retnam_Code/SNDEC_script.R",package = "skewtlmm")
        source(script_path)
    }
    else if (option == "skewN-normal") {
        script_path <- system.file("/Reuben_Retnam_Code/SNnormal_script.R",package = "skewtlmm")
        source(script_path)
        para <- SNLMM.AECM(Y=Y,
                           X=X,
                           Z=Z,
                           N=N,
                           ni=ni,
                           T_i=T_i,
                           be=be_init,
                           la=1,
                           s2=2,
                           Ga=matrix(c(0.35),1,1),
                           tol=1e-4,
                           max.iter=5e5,
                           per=1)
    }
    else if (option == "normal-DEC") {
        script_path <- system.file("/Reuben_Retnam_Code/NLMM_DEC.R",package = "skewtlmm")
        source(script_path)
        para <- NLMM.AECM.DEC(Y=Y,
                              X=X,
                              Z=Z,
                              N=N,
                              ni=ni,
                              T_i=T_i,
                              be=be_init,
                              s2=2,
                              Ga=matrix(c(0.35),1,1),
                              covar_params=c(.5,.9),
                              tol=1e-4,
                              max.iter=5e5,
                              per=1)
    }
    else if (option == "normal-normal") {
        script_path <- system.file("/Reuben_Retnam_Code/NLMM_normal.R",package = "skewtlmm")
        source(script_path)
        para <- NLMM.AECM(Y=Y,
                          X=X,
                          Z=Z,
                          N=N,
                          ni=ni,
                          T_i=T_i,
                          be=be_init,
                          s2=2,
                          Ga=matrix(c(0.35),1,1),
                          tol=1e-4,
                          max.iter=5e5,
                          per=1)
    }
    stopCluster(cl)
    output = para
    return(output)
}
