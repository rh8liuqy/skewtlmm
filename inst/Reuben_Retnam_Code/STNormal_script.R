library(parallel)
library(RcppArmadillo)
library(Rcpp)
#sourceCpp('fastLik_ND.cpp')


calc_expect_1 = function(nu, ni, m, s2, A, sg, sg2, cent, dd, Psi, X, Y, Z) {
  
  s = sqrt(s2)
  tau = pmap_dbl(list(ni, m), ~(nu + .x) / (nu + .y / s2))
  c0 = pmap_dbl(list(A,tau),function(a.i,t.i) a.i * sqrt(t.i) / s)
  c2 = pmap_dbl(list(c0, ni),function(c0.i, ni.i) c0.i * sqrt((nu+ni.i+2)/(nu+ni.i)))
  tmp = pmap_dbl(list(c0, ni), function(c0.i,ni.i) dt(c0.i, df=nu+ni.i) / pt(c0.i, df=nu+ni.i))
  hs80 = pmap_dbl(list(tau, tmp), function(tau.i, tmp.i) sqrt(tau.i) / s * tmp.i)
  hs8 = c0 * tmp
  hs2 = pmap_dbl(list(tau, c2, c0, ni), function(tau.i, c2.i, c0.i, ni.i) tau.i * pt(c2.i, df=nu+ni.i+2) / pt(c0.i, df=nu+ni.i))
  hs3 = sg * A * hs2 + s2 * sg * hs80
  hs4 = sg2 * s2 * (A^2 * hs2 / s2 + 1 + hs8)
  #hs5 = t(t(ub) * hs3) + t(t(vb) * hs2)
  #hs6 = t(t(ub) * hs4) + t(t(vb) * hs3)
  
  #hs7 = hs5 * vb + hs6 * ub + Sb.m * s2
  
  sum.hs2 = sum(hs2)
  sum.hs3 = sum(hs3)
  sum.hs4 = sum(hs4)
  #sum.hs5 = sum(hs5) #rowSums(hs5)
  #sum.hs6 = sum(hs6) #rowSums(hs6)
  #sum.hs7 = sum(hs7) #matrix(rowSums(hs7),q,q) #CHECK 
  XPsiXL = XPsiZL = ZPsiZL = XPsiYL = ZPsiYL = list()
  XPsiXL[[1]] =  hs2[1] * t(X[[1]]) %*% solve(Psi[[1]]) %*% X[[1]]
  XPsiZL[[1]] =   hs3[1] * t(X[[1]]) %*% solve(Psi[[1]]) %*% Z[[1]]
  ZPsiZL[[1]] =   hs4[1] * t(Z[[1]]) %*% solve(Psi[[1]]) %*% Z[[1]]
  XPsiYL[[1]] =   hs2[1] * t(X[[1]]) %*% solve(Psi[[1]]) %*% Y[[1]]
  ZPsiYL[[1]] =   hs3[1] * t(Z[[1]]) %*% solve(Psi[[1]]) %*% Y[[1]]
  Ups1 = pmap(list(hs2, hs3, hs4, cent, dd), function(h2, h3, h4, ce, d) h2*(ce %*% t(ce))+h4 * (d %*% t(d)) - 2*h3*(d %*% t(ce)))
  
  for (i in 2:length(ni)) {
    XPsiXL[[i]] =   hs2[i] * t(X[[i]]) %*% solve(Psi[[i]]) %*% X[[i]]
    XPsiZL[[i]] =    hs3[i] * t(X[[i]]) %*% solve(Psi[[i]]) %*% Z[[i]]
    ZPsiZL[[i]] =    hs4[i] * t(Z[[i]]) %*% solve(Psi[[i]]) %*% Z[[i]]
    XPsiYL[[i]] =    hs2[i] * t(X[[i]]) %*% solve(Psi[[i]]) %*% Y[[i]]
    ZPsiYL[[i]]  =   hs3[i] * t(Z[[i]]) %*% solve(Psi[[i]]) %*% Y[[i]]
    
  }
  XPsiX = Reduce( '+',XPsiXL)
  XPsiZ = Reduce( '+',XPsiZL)
  ZPsiZ = Reduce( '+',ZPsiZL)
  XPsiY = Reduce( '+',XPsiYL)
  ZPsiY = Reduce( '+',ZPsiYL)
  
  #hs2 = hs2, hs3 = hs3, hs4 = hs4, hs5 = hs5, hs6 = hs6, hs7 = hs7, sum.hs2 = sum.hs2, sum.hs3 = sum.hs3, sum.hs4= sum.hs4, sum.hs5= sum.hs5, sum.hs6= sum.hs6, sum.hs7= sum.hs7
  return(list(hs4= hs4,XPsiX = XPsiX,XPsiY = XPsiY, XPsiZ = XPsiZ, ZPsiZ = ZPsiZ, ZPsiY = ZPsiY, Ups1 = Ups1))
}

calc_expect_2 = function(nu, ni, m, s2, A, sg, sg2,ub, vb, Sb.m) {
  
  s = sqrt(s2)
  # tau = (nu + ni) / (nu + m / s2)
  # c0 = A * sqrt(tau) / s
  # c2 = c0 * sqrt((nu+ni+2)/(nu+ni))
  # tmp = dt(c0, df=nu+ni) / pt(c0, df=nu+ni)
  # hs80 = sqrt(tau) / s * tmp
  
  tau = pmap_dbl(list(ni, m), ~(nu + .x) / (nu + .y / s2))
  c0 = pmap_dbl(list(A,tau),function(a.i,t.i) a.i * sqrt(t.i) / s)
  c2 = pmap_dbl(list(c0, ni),function(c0.i, ni.i) c0.i * sqrt((nu+ni.i+2)/(nu+ni.i)))
  tmp = pmap_dbl(list(c0, ni), function(c0.i,ni.i) dt(c0.i, df=nu+ni.i) / pt(c0.i, df=nu+ni.i))
  hs80 = pmap_dbl(list(tau, tmp), function(tau.i, tmp.i) sqrt(tau.i) / s * tmp.i)
  
  hs8 = c0 * tmp
  # hs2 = tau * pt(c2, df=nu+ni+2) / pt(c0, df=nu+ni)
  hs2 = pmap_dbl(list(tau, c2, c0, ni), function(tau.i, c2.i, c0.i, ni.i) tau.i * pt(c2.i, df=nu+ni.i+2) / pt(c0.i, df=nu+ni.i))
  hs3 = pmap_dbl(list(sg, A, hs2, hs80), function(sg.i, A.i, hs2.i, hs80.i) sg.i*A.i*hs2.i + s2*sg.i*hs80.i)#sg * A * hs2 + s2 * sg * hs80
  hs4 = pmap_dbl(list(sg2, A, hs2, hs8), function(sg2.i, A.i, hs2.i, hs8.i) sg2.i *s2* (A.i^2 * hs2.i / s2 + 1+ hs8.i)) #sg2 * s2 * (A^2 * hs2 / s2 + 1 + hs8)
  hs5 = t(t(ub) * hs3) + t(t(vb) * hs2)
  hs6 = t(t(ub) * hs4) + t(t(vb) * hs3)
  
  hs7 = pmap_dbl(list(hs5, vb, hs6, ub, Sb.m), function(hs5.i, vb.i, hs6.i, ub.i, Sb.m.i) hs5.i * vb.i + hs6.i * ub.i + Sb.m.i * s2)#hs5 * vb + hs6 * ub + Sb.m * s2
  
  #sum.hs2 = sum(hs2)
  #sum.hs3 = sum(hs3)
  #sum.hs4 = sum(hs4)
  #sum.hs5 = sum(hs5) #rowSums(hs5)
  #sum.hs6 = sum(hs6) #rowSums(hs6)
  #sum.hs7 = sum(hs7) #matrix(rowSums(hs7),q,q) #CHECK 
  
  
  
  #hs2 = hs2, hs3 = hs3, hs4 = hs4, hs5 = hs5, hs6 = hs6, hs7 = hs7, sum.hs2 = sum.hs2, sum.hs3 = sum.hs3, sum.hs4= sum.hs4, sum.hs5= sum.hs5, sum.hs6= sum.hs6, sum.hs7= sum.hs7
  return(list(hs4 = hs4, hs6 = hs6, hs7 = hs7))
}

STLMM.AECM.ND = function(Y,X,Z,N,ni,T_i, 
                      be, la, s2, Ga, nu, tol = 1e-3, max.iter = 200, per = 100)
{
  #if(length(Y) != N*maxni) return('stop')  
  #if(is.vector(Z) == T) Z = as.matrix(Z)
  begin = proc.time()[1]
  p = ncol(X[[1]]); q = ncol(Z[[1]])
  logli.vector = numeric()
  param.vector = list()
  #X.a = array(X, dim = c(maxni, N, p))
  #Xt.a = aperm(X.a, c(2,1,3))
  #Z.a = array(Z, dim = c(maxni, N, q))
  #Z.uni = matrix(Z.a[,1,], maxni)
  #Y.na = matrix(Y, ncol = maxni, byrow = T)
  
  #----------------------------------------
  diff_vec <<- numeric()
  #diag(maxni)
  #  na.posi = is.na(Y.na)
  #ni = maxni - rowSums(na.posi)
  tgen = lapply(T_i, function(x) abs(outer(x,x,'-')))
  Ip = lapply(T_i, function(x) diag(length(x))) #lapply(T_i, function(t) construct_cov(covar_params, t)) #lapply(ni, function(x) diag(x))
  
  n = sum(ni)
  # ind.na = colSums(2 ^ (1:maxni) * t(na.posi))
  # names(ind.na) = NULL
  # num.cl.na = length(unique(ind.na))
  
  # row.posi = O.list = NULL
  # uni.ind = unique(ind.na)
  # for(i in 1:num.cl.na)
  # {
  #   row.posi[[i]] = which(ind.na == uni.ind[i])
  #   O.list[[i]] = matrix(diag(maxni)[!na.posi[row.posi[[i]][1],],], ncol = maxni)
  # }
  # sub.N = sapply(row.posi, length)
  
  #----------------------------------------
  # local function
  
  sqrt.mt = function(S)
  {
    p = ncol(S)
    if(p == 1) S.sqrt = as.matrix(sqrt(S))
    else
    {
      eig.S = eigen(S)
      S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va),p) %*% t(eig.S$ve)
    }
  }
  
  CML.nu.fn = function(nu, m, A, detL, s2)
  {
    mvt.den = lgamma(.5*(nu+ni))-lgamma(.5*nu)-.5*log(detL)-.5*ni*log(pi*nu*s2)-.5*(nu+ni)*log(1+m/s2/nu)
    Tcdf = log(pt(A*sqrt((nu+ni)/(s2*nu+m)),df=nu+ni))
    val = sum(log(2)+mvt.den+Tcdf)
    return( - val)
  }
  
  CML.covar.fn = function( nu,Ga, Z, FF, delta, cent, s2, tgen, dd)
  {
    # library(purrr)
    # Ip = lapply(T_i, function(x) covar_params[1]^abs(outer(x,x,'-'))^covar_params[2])#lapply(T_i, function(t) construct_cov(covar_params, t))
    # Lambda = pmap(list(Z, Ip), function(z,ip) z %*% Ga %*% t(z) + ip)#Z.uni %*% Ga %*% t(Z.uni) + Ip
    # OLO = Lambda #O %*% Lambda %*% t(O)
    # OLO.inv = lapply(OLO, solve) #solve(OLO)
    # tmp = pmap(list(dd, OLO.inv), function(d, o) t(d) %*% o)#t(Odd)%*% OLO.inv
    # kkk = pmap_dbl(list(tmp, dd), function(t, d) c(1-t %*% d))#c(1-tmp %*% Odd)
    # sg2 = kkk
    # sg = sqrt(kkk) #lapply(kkk, sqrt) #sqrt(kkk)
    # A = pmap_dbl(list(tmp, cent, sg), function(tm, ce, sgi) (tm %*% ce)/sgi)#(tmp %*% ind.cent)/sqrt(kkk)
    # m = pmap_dbl(list(cent, OLO.inv), function(ce, ol.i) t(ce) %*% ol.i %*% ce) #t(cent) %*% OLO %*% cent #colSums((OLO.inv %*% ind.cent) * ind.cent)
    # detL = sapply(OLO, det)#det(OLO)
    
    #A = m  = list()
    A = m  = detL = vector('numeric')
    #mvt.den = Tcdf = 0
    #IV = 0
    val = 0
    for (i in 1:length(Z)) {
      Ip = diag(nrow(Z[[i]])) #covar_params[1]^tgen[[i]]^covar_params[2]
      Lambda = Z[[i]] %*% Ga %*% t(Z[[i]]) + Ip
      OLO.inv = solve(Lambda)
      tmp = t(dd[[i]]) %*% OLO.inv
      kkk = c(1 - tmp %*% dd[[i]])
      A =  (tmp %*% cent[[i]])/sqrt(kkk)
      m = t(cent[[i]]) %*% OLO.inv %*% cent[[i]]
      detL = det(Lambda)
      val = val + log(2) +lgamma(.5*(nu+ni[i]))-lgamma(.5*nu)-.5*log(detL)-.5*ni[i]*log(pi*nu*s2)-.5*(nu+ni[i])*log(1+m/s2/nu) + log(pt(A*sqrt((nu+ni[i])/(s2*nu+m)),df=nu+ni[i]))
    }
    #mvt.den =  lgamma(.5*(nu+ni))-lgamma(.5*nu)-.5*log(detL)-.5*ni*log(pi*nu*s2)-.5*(nu+ni)*log(1+m/s2/nu) 
    #Tcdf = log(pt(A*sqrt((nu+ni)/(s2*nu+m)),df=nu+ni))
    #val = sum(log(2)+mvt.den+Tcdf)
    #val = sum(IV)
    return( - val)
  }
  
  
  #----------------------------------------
  # initial value
  
  #Y.na[na.posi] = 99999
  
  #----------------------------------------
  # observed log-likelihood
  #browser()
  Xb = lapply(X, function(x) x %*% be)#matrix(c(X %*% be), maxni, N)
  cent = pmap(list(Y,Xb), function(y, xb) y - xb)#t(Y.na) - Xb
  Lambda = pmap(list(Y, Z, Ip), function(y,z,ip) z %*% Ga %*% t(z) + ip)#Z.uni %*% Ga %*% t(Z.uni) + Ip
  FF = sqrt.mt(Ga)
  delta = la / sqrt(1 + sum(la^2))
  xi = matrix(c(FF %*% delta),nrow = ncol(Z[[1]])) #as.matrix(c(FF %*% delta))
  dd = map(Z, function(z) z %*% xi) #c(Z.uni %*% xi)
  m = A = detL = sg2 = sg = rep(NA,N)
  
  V = Ga - xi %*% t(xi)
  V.inv = solve(V)
  s = sqrt(s2)
  # for(i in 1:num.cl.na)
  # {
  #index = row.posi[[i]]
  #O = O.list[[i]]
  
  #ub = vb = matrix(NA, q, N)
  #Sb.m = matrix(NA, q^2, N)
  
  #for (i in 1:N) {
  OLO = Lambda #O %*% Lambda %*% t(O)
  OLO.inv = lapply(OLO, solve) #solve(OLO)
  Odd = dd #c(O %*% dd)
  tmp = pmap(list(dd, OLO.inv), function(d, o) t(d) %*% o)#t(Odd)%*% OLO.inv
  ind.cent = cent #O %*% matrix(cent, maxni)
  kkk = pmap_dbl(list(tmp, dd), function(t, d) c(1-t %*% d))#c(1-tmp %*% Odd)
  sg2 = kkk
  sg = sqrt(kkk) #lapply(kkk, sqrt) #sqrt(kkk)
  A = pmap_dbl(list(tmp, cent, sg), function(tm, ce, sgi) (tm %*% ce)/sgi)#(tmp %*% ind.cent)/sqrt(kkk)
  m = pmap_dbl(list(cent, OLO.inv), function(ce, ol.i) t(ce) %*% ol.i %*% ce) #t(cent) %*% OLO %*% cent #colSums((OLO.inv %*% ind.cent) * ind.cent)
  detL = sapply(OLO, det)#det(OLO)
  #OZ = Z.uni #O %*% Z.uni
  OZ = Z
  Sb.inv = pmap_dbl(list(OZ,Ip), function(oz,ip) V.inv + t(oz) %*% solve(ip) %*% oz)   #MAYBE NOT _DBL CHECK THIS 
  Sb = map_dbl(Sb.inv, solve) #solve(Sb.inv)
  Sb.m = c(Sb)
  ub = map_dbl(Sb, function(sb) c(sb %*% V.inv %*% xi))#c(Sb %*% V.inv %*% xi)
  vb = pmap_dbl(list(Sb, OZ, cent, Ip), function(sb, oz, ce, ip) sb %*% t(oz) %*% solve(ip) %*% ce)#Sb %*% t(OZ) %*% ind.cent
  #}
  mvt.den = lgamma(.5*(nu+ni))-lgamma(.5*nu)-.5*log(detL)-.5*ni*log(pi*nu*s2)-.5*(nu+ni)*log(1+m/s2/nu)
  Tcdf = log(pt(A*sqrt((nu+ni)/(s2*nu+m)),df=nu+ni))
  logli.old = sum(log(2)+mvt.den+Tcdf)
  iter = 0
  cat('\n\t\tSTLMM via AECM\n')
  cat(paste(rep('=',60),sep='',collapse=''), '\n')
  cat('iter =', iter, '\tlogli =', logli.old, '\n')
  
  #clusterExport(cl, c('nu', 'idList', 'ni', 'm', 's2', 'A', 'sg', 'sg2', 'cent', 'dd', 'X', 'Y', 'Z', 'calc_expect_1','calc_expect_2'), envir = environment())
  clusterExport(cl, c( 'calc_expect_1','calc_expect_2'), envir = environment())
  repeat
  {
    
    iter = iter + 1
    
    Psi = pmap(list(Z, Ip), function(z, ip) z %*% V %*% t(z) +ip) #Z.uni %*% V %*% t(Z.uni) + Ip
    
    if (iter == 1) {
      
      E1 = list()
      inds = 1:12
      mod_inds = do.call(c, idList[inds])
      
      
      calc_expect_1_ind =  function(i_d) {
        ind_L = idList[[i_d]]
        return(calc_expect_1(nu, ni[ind_L], m[ind_L],  s2, A[ind_L], sg[ind_L], sg2[ind_L], cent[ind_L], dd[ind_L], Psi[ind_L], X[ind_L], Y[ind_L], Z[ind_L]))
      }
      
      estep_upd = foreach(w = inds, .packages = 'purrr') %do% calc_expect_1_ind(w)
      
      
      
      E1$hs4[mod_inds] = do.call(c, map(estep_upd, 'hs4'))
      E1$XPsiZ[inds] = map(estep_upd, 'XPsiZ')
      E1$XPsiY[inds] = map(estep_upd, 'XPsiY')
      E1$XPsiX[inds] = map(estep_upd, 'XPsiX')
      E1$ZPsiZ[inds] = map(estep_upd, 'ZPsiZ')
      E1$ZPsiY[inds] = map(estep_upd, 'ZPsiY')
      E1$Ups1[inds] = map(estep_upd, 'Ups1')
      
      
      
      XPsiZ = Reduce('+',E1$XPsiZ)
      XPsiX = Reduce('+',E1$XPsiX)
      XPsiY = Reduce('+',E1$XPsiY)
      ZPsiZ = Reduce('+',E1$ZPsiZ)
      ZPsiY = Reduce('+',E1$ZPsiY)
      Ups1 = Reduce(c,E1$Ups1)
      
      
    } else {
      
      inds = sample(1:numCores,nwait)
      mod_inds = do.call(c, idList[inds])
      
      
      
      calc_expect_1_ind =  function(i_d) {
        ind_L = idList[[i_d]]
        return(calc_expect_1(nu, ni[ind_L], m[ind_L],  s2, A[ind_L], sg[ind_L], sg2[ind_L], cent[ind_L], dd[ind_L], Psi[ind_L], X[ind_L], Y[ind_L], Z[ind_L]))
      }
      
      estep_upd = foreach(w = inds, .packages = 'purrr') %dopar% calc_expect_1_ind(w)
      
      
      
      E1$hs4[mod_inds] = do.call(c, map(estep_upd, 'hs4'))
      E1$XPsiZ[inds] = map(estep_upd, 'XPsiZ')
      E1$XPsiY[inds] = map(estep_upd, 'XPsiY')
      E1$XPsiX[inds] = map(estep_upd, 'XPsiX')
      E1$ZPsiZ[inds] = map(estep_upd, 'ZPsiZ')
      E1$ZPsiY[inds] = map(estep_upd, 'ZPsiY')
      E1$Ups1[inds] = map(estep_upd, 'Ups1')
      
      
      
      XPsiZ = Reduce('+',E1$XPsiZ)
      XPsiX = Reduce('+',E1$XPsiX)
      XPsiY = Reduce('+',E1$XPsiY)
      ZPsiZ = Reduce('+',E1$ZPsiZ)
      ZPsiY = Reduce('+',E1$ZPsiY)
      Ups1 = Reduce(c,E1$Ups1)
      
    }
    
    # M-step
    
    
    
    #XPsiX = Reduce('+', pmap(list(X, Psi), function( x, ps)  t(x) %*% solve(ps) %*% x))
    #XPsiY= Reduce('+', pmap(list( X, Y, Psi), function( x, y, ps) t(x) %*% solve(ps) %*% y))
    
    
    
    bx.part1 = rbind(cbind(XPsiX, XPsiZ), cbind(t(XPsiZ),ZPsiZ))
    bx.part2 = c(XPsiY, ZPsiY)
    bx = c(solve(bx.part1) %*% bx.part2)
    be = bx[1:p]; xi = bx[p+1:q]
    #browser()
    # Xb = matrix(c(X %*% be), maxni, N)
    # cent = t(Y.na) - Xb
    # dd = c(Z.uni %*% xi)
    Xb = lapply(X, function(x) x %*% be)#matrix(c(X %*% be), maxni, N)•
    cent = pmap(list(Y,Xb), function(y, xb) y - xb)#t(Y.na) - Xb
    dd = map(Z, function(z) z %*% xi) #c(Z.uni %*% xi)
    
    #ppp1 = apply(hs2 * t(cent), 1, rep, times = maxni)
    #ppp2 = matrix(rep(c(cent), each = maxni), maxni^2)
    #Ups1 = Ups2 = ppp1 * ppp2
    #Ups1 = Ups1 + c(dd %*% t(dd)) %*% t(hs4)
    #Ups1 = Ups1 - 2 * apply(hs3 * t(cent), 1, rep, times = maxni) * rep(dd, each = maxni)
    
    ##    Ups1 = with(E1,pmap(list(hs2, hs3, hs4, cent, dd), function(h2, h3, h4, ce, d) h2*(ce %*% t(ce))+h4 * (d %*% t(d)) - 2*h3*(d %*% t(ce))))
    # trPsiUps1 = 0
    # for(i in 1:num.cl.na)
    # {
    #   index = row.posi[[i]]
    #   sub.Ups1 = rowSums(matrix(Ups1[,index],maxni^2))
    #   trPsiUps1 = trPsiUps1 + sum(Psi.oo[[i]] * sub.Ups1)
    # }
    
    s2 = sum(pmap_dbl(list(Psi, Ups1,E1$hs4), function(ps, up,h4) sum(diag(solve(ps) %*% up)) + h4)) / (sum(ni) + N) #c(trPsiUps1) + sum.hs4) / (n + N)
    #remove(Ups1)
    
    #E2 = calc_expect_2(nu, ni, m, s2, A, sg, sg2,ub, vb, Sb.m)
    #browser()
    if (iter == 1) {
      E2 = list()
      inds = 1:12
      mod_inds = do.call(c, idList[inds])
      
      # clusterExport(cl, c('nu', 'm', 's2', 'A', 'sg', 'sg2','ub', 'vb', 'Sb.m'), envir = environment())
      # estep_upd = parLapply(cl, inds, function(i_d) {
      #   ind_L = idList[[i_d]]
      #   return(calc_expect_2(nu, ni[ind_L], m[ind_L], s2, A[ind_L], sg[ind_L], sg2[ind_L],ub[ind_L], vb[ind_L], Sb.m[ind_L]))
      # })
      calc_expect_2_ind = function(i_d) {
        ind_L = idList[[i_d]]
        return(calc_expect_2(nu, ni[ind_L], m[ind_L], s2, A[ind_L], sg[ind_L], sg2[ind_L],ub[ind_L], vb[ind_L], Sb.m[ind_L]))
      }
      estep_upd = foreach(w = inds, .packages = 'purrr') %dopar% calc_expect_2_ind(w)
      
      E2$hs4 = do.call(c, map(estep_upd, 'hs4'))
      E2$hs6 = do.call(c, map(estep_upd, 'hs6'))
      E2$hs7 = do.call(c, map(estep_upd, 'hs7'))
    } else {
      inds = sample(1:numCores,nwait)
      mod_inds = do.call(c, idList[inds])
      
      #clusterExport(cl, c('nu', 'm', 's2', 'A', 'sg', 'sg2','ub', 'vb', 'Sb.m'), envir = environment())
      calc_expect_2_ind = function(i_d) {
        ind_L = idList[[i_d]]
        return(calc_expect_2(nu, ni[ind_L], m[ind_L], s2, A[ind_L], sg[ind_L], sg2[ind_L],ub[ind_L], vb[ind_L], Sb.m[ind_L]))
      }
      
      # estep_upd = parLapply(cl, inds, function(i_d) {
      #   ind_L = idList[[i_d]]
      #   return(calc_expect_2(nu, ni[ind_L], m[ind_L], s2, A[ind_L], sg[ind_L], sg2[ind_L],ub[ind_L], vb[ind_L], Sb.m[ind_L]))
      # })
      estep_upd = foreach(w = inds, .packages = 'purrr') %dopar% calc_expect_2_ind(w)
      
      
      E2$hs4[mod_inds] = do.call(c, map(estep_upd, 'hs4'))
      E2$hs6[mod_inds] = do.call(c, map(estep_upd, 'hs6'))
      E2$hs7[mod_inds] = do.call(c, map(estep_upd, 'hs7'))
    }
    
    ################## E-step
    # tau = (nu + ni) / (nu + m / s2)
    # c0 = A * sqrt(tau) / s
    # c2 = c0 * sqrt((nu+ni+2)/(nu+ni))
    # tmp = dt(c0, df=nu+ni) / pt(c0, df=nu+ni)
    # hs80 = sqrt(tau) / s * tmp
    # hs8 = c0 * tmp
    # hs2 = tau * pt(c2, df=nu+ni+2) / pt(c0, df=nu+ni)
    # hs3 = sg * A * hs2 + s2 * sg * hs80
    # hs4 = sg2 * s2 * (A^2 * hs2 / s2 + 1 + hs8)
    # hs5 = t(t(ub) * hs3) + t(t(vb) * hs2)
    # hs6 = t(t(ub) * hs4) + t(t(vb) * hs3)
    # 
    # hs7 = hs5 * vb + hs6 * ub + Sb.m * s2
    # 
    # sum.hs2 = sum(hs2)
    # sum.hs3 = sum(hs3)
    # sum.hs4 = sum(hs4)
    # sum.hs5 = sum(hs5) #rowSums(hs5)
    # sum.hs6 = sum(hs6) #rowSums(hs6)
    # sum.hs7 = sum(hs7) #matrix(rowSums(hs7),q,q) #CHECK 
    #Ups1 = with(E2,pmap(list(hs2, hs3, hs4, cent, dd), function(h2, h3, h4, ce, d) h2*(ce %*% t(ce))+h4 * (d %*% t(d)) - 2*h3*(d %*% t(ce))))
    ##################
    
    sum.Ups3 = with(E2, sum(E2$hs7) + sum(E2$hs4) * xi %*% t(xi) - sum(E2$hs6) %*% t(xi) - xi %*% t(sum(hs6)))
    
    V = sum.Ups3 / (N * s2)
    
    Ga = V + xi %*% t(xi)
    FF = sqrt.mt(Ga)
    la = c(solve(FF) %*% xi) / c(sqrt(1 - t(xi) %*% solve(Ga) %*% xi))
    
    #Lambda = Z.uni %*% Ga %*% t(Z.uni) + Ip
    Lambda = pmap(list(Y, Z, Ip), function(y,z,ip) z %*% Ga %*% t(z) + ip)#Z.uni %*% Ga %*% t(Z.uni) + Ip
    
    V.inv = solve(V)
    s = sqrt(s2)
    # for(i in 1:num.cl.na)
    # {
    #   index = row.posi[[i]]
    #   O = O.list[[i]]
    #   OLO = O %*% Lambda %*% t(O)
    #   OLO.inv = solve(OLO)
    #   Odd = c(O %*% dd)
    #   tmp = t(Odd)%*%OLO.inv
    #   ind.cent = O %*% matrix(cent[,index], maxni)
    #   kkk = c(1-tmp%*%Odd)
    #   sg2[index] = kkk
    #   sg[index] = sqrt(kkk)
    #   A[index] = (tmp %*% ind.cent)/sqrt(kkk)
    #   m[index] = colSums((OLO.inv %*% ind.cent) * ind.cent)
    #   detL[index] = det(OLO)
    #   OZ = O %*% Z.uni
    #   Sb.inv = V.inv + t(OZ) %*% OZ
    #   Sb = solve(Sb.inv)
    #   Sb.m[,index] = c(Sb)
    #   ub[,index] = c(Sb %*% V.inv %*% xi)
    #   vb[,index] = Sb %*% t(OZ) %*% ind.cent
    # }
    # 
    
    OLO = Lambda #O %*% Lambda %*% t(O)
    OLO.inv = lapply(OLO, solve) #solve(OLO)
    Odd = dd #c(O %*% dd)
    tmp = pmap(list(dd, OLO.inv), function(d, o) t(d) %*% o)#t(Odd)%*% OLO.inv
    ind.cent = cent #O %*% matrix(cent, maxni)
    kkk = pmap_dbl(list(tmp, dd), function(t, d) c(1-t %*% d))#c(1-tmp %*% Odd)
    sg2 = kkk
    sg = sqrt(kkk) #lapply(kkk, sqrt) #sqrt(kkk)
    A = pmap_dbl(list(tmp, cent, sg), function(tm, ce, sgi) (tm %*% ce)/sgi)#(tmp %*% ind.cent)/sqrt(kkk)
    m = pmap_dbl(list(cent, OLO.inv), function(ce, ol.i) t(ce) %*% ol.i %*% ce) #t(cent) %*% OLO %*% cent #colSums((OLO.inv %*% ind.cent) * ind.cent)
    detL = sapply(OLO, det)#det(OLO)
    #OZ = Z.uni #O %*% Z.uni
    OZ = Z
    Sb.inv = pmap_dbl(list(OZ,Ip), function(oz,ip) V.inv + t(oz) %*% solve(ip) %*% oz)   #MAYBE NOT _DBL CHECK THIS 
    Sb = map_dbl(Sb.inv, solve) #solve(Sb.inv)
    Sb.m = c(Sb)
    ub = map_dbl(Sb, function(sb) c(sb %*% V.inv %*% xi))#c(Sb %*% V.inv %*% xi)
    vb = pmap_dbl(list(Sb, OZ, cent, Ip), function(sb, oz, ce, ip) sb %*% t(oz) %*% solve(ip) %*% ce)#Sb %*% t(OZ) %*% ind.cent
    
    nu.optim = optim(par = nu, fn = CML.nu.fn, method = "L-BFGS-B", lower = 4, upper = 15, m=m, A=A, detL=detL, s2=s2, control = list(factr = 1e16))
    nu = nu.optim$par
    
    #if ((iter %% 5) == 0 | iter == 1) {
      #rho.optim = optim(par = covar_params, fn = fastLikelihood, method = 'L-BFGS-B', lower = c(0.01, 0.2), upper = c(0.99,2), nu = nu,ni = ni, Ga = Ga, Z =Z, cent = cent, delta = delta, s2 = s2, T_i = tgen, dd = dd, N = length(ni), control = list(factr = 1e16))
      #covar_params = c(1,0)#rho.optim$par
      #logli.new = -rho.optim$value
    #} else {
      logli.new = -CML.covar.fn(nu,Ga, Z, FF, delta, cent, s2, tgen, dd)#-fastLikelihood( nu = nu,ni = ni, Ga = Ga, Z =Z, cent = cent, delta = delta, s2 = s2, T_i = tgen, dd = dd, N = length(ni))
    #}
    #Ip = lapply(T_i, function(x) covar_params[1]^abs(outer(x,x,'-'))^covar_params[2]) #lapply(T_i, function(t) construct_cov(covar_params, t))
    
    
    
    
    diff = logli.new - logli.old
    diff_vec <<- c(diff_vec, diff)
    if(iter %% per == 0) cat('iter =', iter, '\tlogli =', logli.new, '\tdiff =', diff, '\n')
    if(abs(diff) < tol | iter == max.iter) break
    logli.old = logli.new
    logli.vector = c(logli.vector, logli.new)
    param.vector[[iter]] = c(be, delta, xi, s2, Ga, FF, nu)
  }
  num.para = p+q+1+q*(q+1)/2+1
  aic = -2 * logli.new + 2      * num.para
  bic = -2 * logli.new + log(N) * num.para
  
  model.inf = c(num.para, logli.new, aic, bic)
  names(model.inf) = c('num.para', 'logli', 'AIC', 'BIC')
  delta = c(solve(FF) %*% xi)
  para = list(model.inf = model.inf, beta = be, lambda = cbind(la,delta,xi), s2 = c(s=s,s2=s2), Gamma = Ga, F = FF, nu = nu, time = c(proc.time()[1]-begin, iter, logli.new, logli.vector), param.vector)
  cat(paste(rep('=',60),sep='',collapse=''), '\n')
  cat('Converged observed logli:', logli.new, ';\t Iterations:', iter, '\n')
  #  end = proc.time()[1]
  #  cat("It tooks", end - begin, "seconds\n\n")
  
  #---------------------------------------------------
  # Estimation of random effects (Method: E-step)
  
  # tau = (nu + ni) / (nu + m / s2)
  # c0 = A * sqrt(tau) / s
  # c2 = c0 * sqrt((nu+ni+2)/(nu+ni))
  # cm2 = c0 * sqrt((nu+ni-2)/(nu+ni))
  # 
  # coef = sg * A * (1 + dt(cm2, df=nu+ni-2) / pt(c0, df=nu+ni) / cm2)
  # ran.eff = t(t(ub) * coef) + vb
  
  #---------------------------------------------------
  #fitted values and residuals
  
  # y.fit = Xb + Z.uni %*% ran.eff
  # resi = Y - y.fit
  
  #---------------------------------------------------
  # prediction
  # 
  # pred = t(Y.na)
  # for(i in 1:num.cl.na)
  # {
  #   index = row.posi[[i]]
  #   ind.cent = matrix(cent[,index],maxni)
  #   ind.Xb = matrix(Xb[,index],maxni)
  #   O = O.list[[i]]
  #   OO = t(O) %*% O
  #   w = ind.Xb + OO %*% ind.cent
  #   W = (diag(maxni) - OO) %*% Z.uni
  #   pred[,index] = w + W %*% ran.eff[,index]
  # }
  
  #---------------------------------------------------
  
  outl = cbind(A=A, Delta=m, ni=ni)
  end = proc.time()[1]
  cat('It took', end - begin, 'seconds\n')
  cat('MLE:\n')
  print(para)
  return(list(para))
}
