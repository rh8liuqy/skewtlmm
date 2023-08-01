
NLMM.AECM = function(Y,X,Z,N,ni,T_i, 
                      be, s2, Ga, tol = 1e-3, max.iter = 200, per = 100)
{
  #if(length(Y) != N*maxni) return('stop')  
  #if(is.vector(Z) == T) Z = as.matrix(Z)
  la = 0
  begin = proc.time()[1]
  p = ncol(X[[1]]); q = ncol(Z[[1]])
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
  Ip = lapply(ni, function(x) diag(x))
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
  
  CML.fn = function( m, A, detL, s2)
  {
    mvn.den = -.5*log(detL)-.5*ni*log(2*pi*s2)-.5*m/s2
    Ncdf = log(pnorm(A/sqrt(s2))+0.001)
    logli.new = sum(log(2)+mvn.den+Ncdf)
    return( -logli.new)
  }
  
  CML.covar.fn = function(covar_params,Ga, Z, FF, delta, s2, T_i)
  {
    Ip = lapply(T_i, function(t) construct_cov(covar_params, t))
    xi = matrix(c(FF %*% delta),nrow = ncol(Z[[1]])) #as.matrix(c(FF %*% delta))
    dd = map(Z, function(z) z %*% xi) #c(Z.uni %*% xi)
    Lambda = pmap(list(Z, Ip), function(z,ip) z %*% Ga %*% t(z) + ip)#Z.uni %*% Ga %*% t(Z.uni) + Ip
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
    
    
    mvn.den = -.5*log(detL)-.5*ni*log(2*pi*s2)-.5*m/s2
    Ncdf = log(pnorm(A/sqrt(s2)))
    logli.new = sum(log(2)+mvn.den+Ncdf)
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
  Sb.inv = map_dbl(OZ, function(oz) V.inv + t(oz) %*% oz)   #MAYBE NOT _DBL CHECK THIS 
  Sb = map_dbl(Sb.inv, solve) #solve(Sb.inv)
  Sb.m = c(Sb)
  ub = map_dbl(Sb, function(sb) c(sb %*% V.inv %*% xi))#c(Sb %*% V.inv %*% xi)
  vb = pmap_dbl(list(Sb, OZ, cent), function(sb, oz, ce) sb %*% t(oz) %*% ce)#Sb %*% t(OZ) %*% ind.cent
  #}
  #mvt.den = lgamma(.5*(nu+ni))-lgamma(.5*nu)-.5*log(detL)-.5*ni*log(pi*nu*s2)-.5*(nu+ni)*log(1+m/s2/nu)
  #Tcdf = log(pt(A*sqrt((nu+ni)/(s2*nu+m)),df=nu+ni))
  logli.old = -CML.fn( m, A, detL, s2)
  #logli.old = sum(log(2)+mvt.den+Tcdf)
  iter = 0
  
  cat('\n\t\tSTLMM via AECM\n')
  cat(paste(rep('=',60),sep='',collapse=''), '\n')
  cat('iter =', iter, '\tlogli =', logli.old, '\n')
  
  repeat
  {
    #browser()
    iter = iter + 1
    # tau = (nu + ni) / (nu + m / s2)
    # c0 = A * sqrt(tau) / s
    # c2 = c0 * sqrt((nu+ni+2)/(nu+ni))
    # 
    # # E-step
    # 
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
    hs6 = (dnorm(A/s)+0.001) / (pnorm(A/s)+0.001)
    hs1 = sg * A + s * sg * hs6
    hs2 = sg2 * s2 + sg * A * hs1
    hs3 = t(t(ub) * hs1) + vb
    #browser()
    hs4 = matrix(ub * hs2 + vb * hs1,ncol = 1)
    hs5 = vb * hs3 + ub * hs4 + Sb.m * s2 #apply(hs3, 2, rep, times = q) * matrix(rep(c(vb), each = q), q^2) +
    #apply(hs4, 2, rep, times = q) * matrix(rep(c(ub), each = q), q^2) + Sb.m * s2
    
    #apply(hs5, 2, rep, times = q) * matrix(rep(c(vb), each = q), q^2) +
    #apply(hs6, 2, rep, times = q) * matrix(rep(c(ub), each = q), q^2) + Sb.m * s2
    
    sum.hs2 = sum(hs2)
    sum.hs3 = sum(hs3)
    sum.hs4 = sum(hs4)
    sum.hs5 = sum(hs5) #rowSums(hs5)
    #sum.hs6 = sum(hs6) #rowSums(hs6)
    #sum.hs7 = sum(hs7) #matrix(rowSums(hs7),q,q) #CHECK 
    
    # M-step
    
    Psi = pmap(list(Z, Ip), function(z, ip) z %*% V %*% t(z) +ip) #Z.uni %*% V %*% t(Z.uni) + Ip
    #XPsiX = XPsiZ = ZPsiZ = XPsiY = ZPsiY = 0
    #Psi.oo = as.list(numeric(num.cl.na))
    # for(i in 1:num.cl.na)
    # {
    #   index = row.posi[[i]]
    #   O = O.list[[i]]
    #   Psi.oo[[i]] = t(O) %*% solve(O %*% Psi %*% t(O)) %*% O
    # for(j in index)
    # {
    #   XPsiX = XPsiX + hs2[j] * t(X.a[,j,]) %*% Psi.oo[[i]] %*% X.a[,j,]
    #   XPsiZ = XPsiZ + hs3[j] * t(X.a[,j,]) %*% Psi.oo[[i]] %*% Z.uni
    #   ZPsiZ = ZPsiZ + hs4[j] * t(Z.uni) %*% Psi.oo[[i]] %*% Z.uni
    #   XPsiY = XPsiY + hs2[j] * t(X.a[,j,]) %*% Psi.oo[[i]] %*% Y.na[j,]
    #   ZPsiY = ZPsiY + hs3[j] * t(Z.uni) %*% Psi.oo[[i]] %*% Y.na[j,]
    # }
    
    XPsiX = Reduce('+', pmap(list(X, Psi), function( x, ps)  t(x) %*% solve(ps) %*% x))
    XPsiZ = Reduce('+', pmap(list(hs1, X, Psi,Z), function(h1, x, ps,z) h1 * t(x) %*% solve(ps) %*% z))
    ZPsiZ = Reduce('+', pmap(list(hs2, Z, Psi), function(h2, z, ps) h2 * t(z) %*% solve(ps) %*% z))
    XPsiY= Reduce('+', pmap(list( X, Y, Psi), function( x, y, ps) t(x) %*% solve(ps) %*% y))
    ZPsiY = Reduce('+', pmap(list(hs1, Z, Y, Psi), function(h1, z, y, ps) h1 * t(z) %*% solve(ps) %*% y))
    #}
    
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
    Ups1 = pmap(list(hs1, hs2, cent, dd), function(h1, h2, ce, d) (ce %*% t(ce))+h2*(d %*% t(d)) - 2*h1*(d %*% t(ce)))
    # trPsiUps1 = 0
    # for(i in 1:num.cl.na)
    # {
    #   index = row.posi[[i]]
    #   sub.Ups1 = rowSums(matrix(Ups1[,index],maxni^2))
    #   trPsiUps1 = trPsiUps1 + sum(Psi.oo[[i]] * sub.Ups1)
    # }
    
    s2 = sum(pmap_dbl(list(Psi, Ups1,hs2), function(ps, up,h2) sum(diag(solve(ps) %*% up)) + h2)) / (sum(ni) + N) #c(trPsiUps1) + sum.hs4) / (n + N)
    #remove(Ups1)
    
    #sum.Ups3 = sum.hs7 + sum.hs4 * xi %*% t(xi) - sum.hs6 %*% t(xi) - xi %*% t(sum.hs6)
    sum.Ups3 = sum.hs5 + sum.hs2 * xi %*% t(xi) - sum.hs4 %*% t(xi) - xi %*% t(sum.hs4)
    V = sum.Ups3 / (N * s2)
    
    Ga = V + xi %*% t(xi)
    FF = sqrt.mt(Ga)
    la = 0#c(solve(FF) %*% xi) / c(sqrt(1 - t(xi) %*% solve(Ga) %*% xi))
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
    Sb.inv = map_dbl(OZ, function(oz) V.inv + t(oz) %*% oz)   #MAYBE NOT _DBL CHECK THIS 
    Sb = map_dbl(Sb.inv, solve) #solve(Sb.inv)
    Sb.m = c(Sb)
    ub = map_dbl(Sb, function(sb) c(sb %*% V.inv %*% xi))#c(Sb %*% V.inv %*% xi)
    vb = pmap_dbl(list(Sb, OZ, cent), function(sb, oz, ce) sb %*% t(oz) %*% ce)#Sb %*% t(OZ) %*% ind.cent
    
    #nu.optim = optim(par = nu, fn = CML.nu.fn, method = "L-BFGS-B", lower = 4, upper = 20, m=m, A=A, detL=detL, s2=s2)
    #nu = nu.optim$par
    
    
    
    logli.new = -CML.fn( m, A, detL, s2)#-nu.optim$value
    diff = logli.new - logli.old
    diff_vec <<- c(diff_vec, diff)
    if(iter %% per == 0) cat('iter =', iter, '\tlogli =', logli.new, '\tdiff =', diff, '\n')
    if((diff) < tol | iter == max.iter) break
    logli.old = logli.new
  }
  num.para = p+q+1+q*(q+1)/2+1
  aic = -2 * logli.new + 2      * num.para
  bic = -2 * logli.new + log(N) * num.para
  model.inf = c(num.para, logli.new, aic, bic)
  names(model.inf) = c('num.para', 'logli', 'AIC', 'BIC')
  delta = c(solve(FF) %*% xi)
  para = list(model.inf = model.inf, beta = be, lambda = cbind(la,delta,xi), s2 = c(s=s,s2=s2), Gamma = Ga, F = FF)
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
