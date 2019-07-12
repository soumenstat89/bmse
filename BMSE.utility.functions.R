
required.packages = c('mvtnorm', 'MCMCpack', 'abind', 'ggplot2', 'matrixStats',
                     'Brobdingnag', 'parallel', 'snowfall', 'rlecuyer', 'stargazer')
notinstalled = setdiff(required.packages,  installed.packages()[,1])
if (length(notinstalled) > 0) install.packages(notinstalled, dependencies = T)
sapply(required.packages, require, character.only = TRUE)

#====================================================================
# define utility functions
#====================================================================
Mode <- function(x) 
{
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

logit = function(x) { log(x) - log(1 - x) }
expit = function(x) { ifelse(x < 0, exp(x) /(1 + exp(x)),  1/(1 + exp(-x))) }
lowerQuantile = function(z) {  quantile(z, probs = 0.025)  }
upperQuantile = function(z) {  quantile(z, probs = 0.975)  }
medQuantile = function(z) {  quantile(z, probs = 0.5)  }
bias = function(x){  mean((x[-1]-x[1]), na.rm =T)}
relbiaspercentage = function(x){  100*(mean(x[-1], na.rm =T)-x[1])/x[1]}
rmse = function(x){  sqrt(mean((x[-1]-x[1])^2, na.rm =T))}
sim.stats = function(x, hpdprob = 0.95)
{
  bias =   mean(x[-1], na.rm =T) - x[1]
  relbiaspercentage = 100*(mean(x[-1], na.rm =T)-x[1])/x[1] 
  rmse =  sqrt(mean((x[-1]-x[1])^2, na.rm =T))
  hpd = HPDinterval(as.mcmc(x[-1]), prob = hpdprob)
  return(c(bias, relbiaspercentage, rmse, hpd))
}

e2dist = function(x,y) # m x n # 'x' is matrix mx2, 'y' is matrix nx2
{
  ind = sort(rep(1:nrow(y), nrow(x))) 
  dvec = sqrt((x[,1] - y[ind,1]) ^ 2 + (x[,2] - y[ind,2]) ^ 2) 
  dmat = matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F) 
  return(dmat)
}

sumvec = function(a){ s = a[1]; for(i in 2:length(a)) s = s + a[i]; return(s)}
scorefn = function(x, mu, Sigma){ t(x-mu)%*% solve(Sigma) %*% (x-mu)}
dmvtruncnorm = function(x, mu, Sigma, conf.lev, logg) { dmvnorm(x, mean = mu, sigma = Sigma, log = logg) - log(conf.lev)}

logMeanExp = function(x, na.rm = T)
{
  sumna = sum(is.na(x))
  if(sumna == 0){  a = logSumExp(x) - log(length(x))}
  if(sumna > 0){  a = logSumExp(x, na.rm = T) - log(length(x) - sumna)}
  return(a)
}  

gdmean = function(logh.chain)
{
  C = mean(logh.chain, na.rm = T)
  est = - (C + logMeanExp(logh.chain - C))
  return(est)
}

logpimatfn = function(logitbasepar, s.left, traploc, sex, logsigmam, logsigmaf, Msexsigma)
{
  if(Msexsigma == 1){  logsigma = rep(logsigmaf, nrow(s.left)); logsigma[sex == 1] = logsigmam }
  if(Msexsigma == 0){  logsigma = logsigmam }
  D1 = e2dist(s.left, traploc) # M x ntrap
  logpimat = log(expit(logitbasepar)) - D1 * D1 / (2 * (exp(logsigma) ^ 2)) # Division is always column wise
  return(logpimat)
}

logpmatfn = function(logitp0, s.left, traploc, sex, logsigmam, logsigmaf, Msexsigma)
{
  if(Msexsigma == 1){  logsigma = rep(logsigmaf, nrow(s.left));  logsigma[sex == 1] = logsigmam }
  if(Msexsigma == 0){  logsigma = logsigmam }
  D1 = e2dist(s.left, traploc)
  logpmat = log(expit(logitp0)) - D1 * D1 / (2 * (exp(logsigma) ^ 2)) # Division is always column wise
  return(logpmat)
}

waicfn = function(loglik.mat) # dim(loglik.mat) = (ndraws - burnin) x M 
{
  t1 = sum(apply(loglik.mat, 2, logMeanExp, na.rm = T)) # scalar
  u1 = apply(loglik.mat, 2, mean, na.rm = T) # M x 1
  u2 = abs(t(t(loglik.mat) - u1)) # (ndraws - burnin) x M 
  t2 = sum(u1) # scalar
  pd1 =  -2 * t2 + 2 * t1 # scalar
  pd2 =  sum(apply(loglik.mat, 2, var, na.rm = T))  # scalar
  pd3 =  sum(apply(u2, 2, mean, na.rm = T)) # scalar 
  waic1 = -2 * t1 + 2 * pd1
  waic2 = -2 * t1 + 2 * pd2
  waic3 = -2 * t1 + 2 * pd3
  return(c(waic1, waic2, waic3, pd1, pd2, pd3))
}

dicfn = function(loglik, loglik.pt) # (ndraws - burnin) x 1 
{
  Dhat = -2 * loglik.pt
  Dbar = -2 * mean(loglik, nar.rm =T)
  pd1 =  Dbar - Dhat # 2 * (loglik.pt - mean(loglik, nar.rm =T))
  pd2 = 2 * var(loglik)
  dic1 = Dhat + 2 * pd1 #  -2 * loglik.pt + 2 * pd1
  dic2 = Dhat + 2 * pd2 #  -2 * loglik.pt + 2 * pd2
  return(c(dic1, dic2, pd1, pd2))
} 
#==============================================
# Dey et al. (2019) Partial ID - log - likelihood function
#==============================================

logLfn.piscr = function(logitphi, logpimat, z, yi00, n0mat, n0vec, J)
{
  phi = expit(logitphi)
  A1 =  log(phi ^ yi00) + log((1 - phi) ^ (2*n0vec - yi00)) # M x 1
  pimat = exp(logpimat)
  pimat[pimat < .Machine$double.xmin] = .Machine$double.xmin
  A2 = rowSums(log(pimat^n0mat), na.rm = T) # M x 1
  A3 = rowSums(log((1 - phi* (2 - phi) * pimat)^(J-n0mat)), na.rm = T) # M x 1
  bignum = -.Machine$double.min.exp; A1[A1 == -Inf] = -bignum; A2[A2 == -Inf] = -bignum; A3[A3 == -Inf] = -bignum;  
  return(sum(z*(A1 + A2 + A3)))
}

logLfn.piscr.i.waic = function(logitphi, logpimat, z, yi00, n0mat, n0vec, J)
{
  phi = expit(logitphi)
  A1 =  log(phi ^ yi00) + log((1 - phi) ^ (2*n0vec - yi00)) # M x 1
  pimat = exp(logpimat)
  pimat[pimat < .Machine$double.xmin] = .Machine$double.xmin
  A2 = rowSums(log(pimat^n0mat), na.rm = T) # M x 1
  A3 = rowSums(log((1 - phi* (2 - phi) * pimat)^(J-n0mat)), na.rm = T) # M x 1
  bignum = -.Machine$double.min.exp; A1[A1 == -Inf] = -bignum; A2[A2 == -Inf] = -bignum; A3[A3 == -Inf] = -bignum;
  return(z*(A1 + A2 + A3))
}

logLfn.piscr.u0zs = function(par, data1)
{
  psi = par$psi; phi = par$phi; omega0 = par$omega0 
  Msexsigma = par$Msexsigma;
  if(Msexsigma == 1)
  {
    sigmam = par$sigmam; 
    sigmaf = par$sigmaf; 
  }
  if(Msexsigma == 0) {sigma = par$sigma} 
  L = unlist(par$L); 
  
  left = data1$left; yli00 = rowSums(left); right = data1$right; yri00 = rowSums(right);
  right.star = right[order(L),,]; yri00.star = yri00[order(L)]; yi00 = yli00 + yri00.star # giving total no. of captures for each captured individual
  n0mat = apply(left + right.star, c(1,2), function(a){sum(a > 0)}) # trapwise count of recorded entries for each of the M individual
  n0vec = rowSums(n0mat) # giving total count of recorded entries for each of the M individual
  zero.guys = yi00 == 0 ; num0 = sum(zero.guys)# M x 1 logical vector indicating captured or not
  numcap = sum(yi00 > 0); 
  
  numl = data1$numl; traploc = data1$traploc; 
  sgridloc = data1$sgridloc; nG = nrow(sgridloc)
  ntrap = nrow(traploc); J = data1$J; M = data1$M;
  bignum = -.Machine$double.min.exp; KnG = ntrap *nG
  
  D.t = e2dist(traploc, sgridloc) #  ntrap x nG # Transpose of the usal distance matrix
  AA =  log(phi^yi00) + log((1 - phi)^(2*n0vec - yi00)) # M x 1
  #======================================= 
  if(Msexsigma == 1)
  {
    theta = par$theta; 
    sex.obs = data1$sex.obs # numl x 1
    logpimatm.t  = log(omega0) - D.t * D.t / (2 * (sigmam  ^ 2)) # ntrap x nG # Transpose of the usual pimatm
    logpimatf.t  = log(omega0) - D.t * D.t / (2 * (sigmaf  ^ 2)) # ntrap x nG # Transpose of the usual pimatf
    pimatm.t = exp(logpimatm.t)
    pimatm.t[pimatm.t < .Machine$double.xmin] = .Machine$double.xmin
    pimatf.t = exp(logpimatf.t)
    pimatf.t[pimatf.t < .Machine$double.xmin] = .Machine$double.xmin
    U1 = U2 = U3 = U4 = matrix(NA, nrow = ntrap, ncol = nG) # ntrap x nG
    
    #--------------------------------------------------
    # For the captured and recorded gender individuals
    #--------------------------------------------------
    
    A1 = A2 = matrix(NA, numl, nG) # numl x nG
    for(ind in 1:numl)
    {
      if(sex.obs[ind] == 1)
      {
        U1[1:KnG] = log(pimatm.t[1:KnG]^rep(n0mat[ind,], nG))
        U2[1:KnG] = log((1 - phi* (2 - phi) * pimatm.t)[1:KnG] ^ rep(J-n0mat[ind,], nG))
      }
      if(sex.obs[ind] == 0)
      {
        U1[1:KnG] = log(pimatf.t[1:KnG]^rep(n0mat[ind,], nG))
        U2[1:KnG] = log((1 - phi* (2 - phi) * pimatf.t)[1:KnG] ^ rep(J-n0mat[ind,], nG))
      }
      A1[ind,] = colSums(U1, na.rm = T) # nG x 1 # Summing over the traps.
      A2[ind,] = colSums(U2, na.rm = T) # nG x 1 # Summing over the traps.
    }
    
    A1[A1 == -Inf] = -bignum; A2[A2 == -Inf] = -bignum; 
    
    part1 = log(psi) + log(theta^sex.obs) + log((1-theta)^(1-sex.obs)) +
      AA[1:numl] + apply(A1+A2, 1, logMeanExp, na.rm = T) # numl x 1
    
    #--------------------------------------------------
    # For the captured but missing gender individuals
    #--------------------------------------------------
    if(numcap > numl){
      III = c(1:M)[yi00 > 0 & c(1:M) > numl] # (numcap - numl) x 1
      B1 = B2 = B3 = B4 = matrix(NA, length(III), nG) # (numcap - numl) x nG
      i = 0
      for(ind in III)
      {
        i = i + 1
        U1[1:KnG] = log(pimatm.t[1:KnG]^rep(n0mat[ind,], nG))
        U2[1:KnG] = log((1 - phi* (2 - phi) * pimatm.t)[1:KnG] ^ rep(J-n0mat[ind,], nG))
        U3[1:KnG] = log(pimatf.t[1:KnG]^rep(n0mat[ind,], nG)) 
        U4[1:KnG] = log((1 - phi* (2 - phi) * pimatf.t)[1:KnG] ^ rep(J-n0mat[ind,], nG))
        
        B1[i,] = colSums(U1, na.rm = T)  # Summing over the traps.
        B2[i,] = colSums(U2, na.rm = T)  # Summing over the traps.
        B3[i,] = colSums(U3, na.rm = T)  # Summing over the traps.
        B4[i,] = colSums(U4, na.rm = T)  # Summing over the traps.
      }
      
      B1[B1 == -Inf] = -bignum; B2[B2 == -Inf] = -bignum; 
      B3[B3 == -Inf] = -bignum; B4[B4 == -Inf] = -bignum; 
   
      part2 = log(psi) + AA[III] +  log(theta * exp(apply(B1+B2, 1, logMeanExp, na.rm = T))+ 
                                        (1-theta) * exp(apply(B3+B4, 1, logMeanExp, na.rm = T))) # (numcap - numl) x 1
                                
    } # end of if(numcap > numl)
    #----------------------------------------------------
    # For the uncaptured and missing gender individuals
    #----------------------------------------------------
    U1[1:KnG] = log((1 - phi* (2 - phi) * pimatm.t)[1:KnG] ^ rep(J, KnG))
    U2[1:KnG] = log((1 - phi* (2 - phi) * pimatf.t)[1:KnG] ^ rep(J, KnG))
    
    C1 = colSums(U1, na.rm = T) # nG x 1 # Summing over the traps.
    C2 = colSums(U2, na.rm = T) # nG x 1 # Summing over the traps.
    C1[C1 == -Inf] = -bignum; C2[C2 == -Inf] = -bignum
    
    part3 = log(psi*(theta*mean(exp(C1),na.rm=T)+(1-theta)*mean(exp(C2),na.rm=T)) + (1-psi)) # scalar
    
    #----------------------------------
    if(numcap > numl){loglik = sum(part1) + sum(part2) + (M-numcap) * part3}
    if(numcap == numl){loglik = sum(part1) + (M-numcap) * part3}
    return(loglik) 
  } # end of  if(Msexsigma == 1)
  #=============================== 
  if(Msexsigma == 0)
  { 
    logpimat.t = log(omega0) - D.t * D.t / (2 * (sigma  ^ 2)) # ntrap x nG # Transpose of the usual pimat
    pimat.t = exp(logpimat.t)
    pimat.t[pimat.t < .Machine$double.xmin] = .Machine$double.xmin
    #-------------------------------
    # For the captured individuals
    #-------------------------------
    III = c(1:M)[yi00 > 0] # numcap x 1
    
    B1 = B2 = matrix(NA, numcap, nG) # numcap x nG
    U1 = U2 = matrix(NA, ntrap, nG) # ntrap x nG
    i = 0
    for(ind in III)
    {
      i = i+1
      U1[1:KnG] = log(pimat.t[1:KnG]^rep(n0mat[ind,], nG))
      U2[1:KnG] = log((1 - phi* (2 - phi) * pimat.t)[1:KnG] ^ rep(J-n0mat[ind,], nG))
      B1[i,] = colSums(U1, na.rm = T) # nG x 1 # Summing over the traps.
      B2[i,] = colSums(U2, na.rm = T) # nG x 1 # Summing over the traps.
    }
    B1[B1 == -Inf] = -bignum; B2[B2 == -Inf] = -bignum; 
    part1 = log(psi) + AA[III] + apply(B1 + B2, 1, logMeanExp, na.rm = T) # numcap x 1
    
    #---------------------------------
    # For the uncaptured individuals
    #---------------------------------
    U1[1:KnG] = log((1 - phi* (2 - phi) * pimat.t)[1:KnG] ^ rep(J, KnG))
    
    C1 = colSums(U1, na.rm = T) # nG x 1 # Summing over the traps.
    C1[C1 == -Inf] = -bignum; 
    part2 = log(psi * mean(exp(C1), na.rm = T) + (1-psi)) # scalar
    #----------------------------------     
    loglik = sum(part1) + (M-numcap) * part2
    return(loglik) 
  } # end of  if(Msexsigma == 0)
  
}
#==============================================
# Royle (2015) Partial ID: log-likelihood function
#==============================================
logLfn.rpiscr = function(logpmat, left2d, right2d, z, J) # Used to update L
{ 
  part1 = rowSums(dbinom(left2d, J, exp(logpmat), log = T)) # M x 1
  part2 = rowSums(dbinom(right2d, J, exp(logpmat), log = T)) # M x 1
  bignum = -.Machine$double.min.exp; part1[part1 == -Inf] = -bignum; part2[part2 == -Inf] = -bignum; 
  loglik = sum(z*(part1 + part2))
  return(loglik)
}

logLfn.rpiscr.dic = function(logpmat, left2d, right2d, z, J) 
{ 
  caphist2d = left2d + right2d
  pmat = exp(logpmat)
  pmat[pmat < .Machine$double.xmin] = .Machine$double.xmin
  part1 = rowSums(log(pmat^caphist2d)) # M x 1 
  part2 = rowSums(log((1 - pmat) ^ (2*J - caphist2d))) # M x 1 
  bignum = -.Machine$double.min.exp; part1[part1 == -Inf] = -bignum; part2[part2 == -Inf] = -bignum; 
  loglik = sum(z*(part1 + part2))
  return(loglik)
}

logLfn.rpiscr.i.waic = function(logpmat, left2d, right2d, z, J) 
{
  caphist2d = left2d + right2d
  pmat = exp(logpmat)
  pmat[pmat < .Machine$double.xmin] = .Machine$double.xmin
  part1 = rowSums(log(pmat^caphist2d)) # M x 1 
  part2 = rowSums(log((1 - pmat) ^ (2*J - caphist2d))) # M x 1 
  bignum = -.Machine$double.min.exp; part1[part1 == -Inf] = -bignum; part2[part2 == -Inf] = -bignum;  
  loglik = z*(part1 + part2)
  return(loglik)
}

logLfn.rpiscr.u0zs = function(par, data1)
{
  psi = par$psi; p0 = par$p0 
  Msexsigma = par$Msexsigma;
  if(Msexsigma == 1)
  {
    sigmam = par$sigmam; 
    sigmaf = par$sigmaf; 
  }
  if(Msexsigma == 0) {sigma = par$sigma} 
  L = unlist(par$L); 
  
  left = data1$left; right = data1$right
  left2d = data1$left2d; right2d = data1$right2d
  right2d.star = right2d[order(L),];
  
  caphist2d = left2d + right2d.star
  yi00 = rowSums(left2d) + rowSums(right2d.star)
  zero.guys = yi00 == 0 ; num0 = sum(zero.guys) # M x 1 logical vector indicating captured or not
  numcap = sum(yi00 > 0)
  
  numl = data1$numl; traploc = data1$traploc; 
  sgridloc = data1$sgridloc; nG = nrow(sgridloc)
  ntrap = nrow(traploc); J = data1$J; M = data1$M;
  bignum = -.Machine$double.min.exp; KnG = ntrap*nG
  
  D.t = e2dist(traploc, sgridloc) #  ntrap x nG # Transpose of the usal distance matrix
  
  if(Msexsigma == 1)
  {
    theta = par$theta 
    sex.obs = data1$sex.obs # numl x 1
    
    logpmatm.t  = log(p0) - D.t * D.t / (2 * (sigmam  ^ 2)) # ntrap x nG # Transpose of the usual pmatm
    logpmatf.t  = log(p0) - D.t * D.t / (2 * (sigmaf  ^ 2)) # ntrap x nG # Transpose of the usual pmatf
    pmatm.t = exp(logpmatm.t)
    pmatm.t[pmatm.t < .Machine$double.xmin] = .Machine$double.xmin
    pmatf.t = exp(logpmatf.t)
    pmatf.t[pmatf.t < .Machine$double.xmin] = .Machine$double.xmin
    
    U1 = U2 = U3 = U4 = matrix(NA, nrow = ntrap, ncol = nG) # ntrap x nG
   
    #--------------------------------------------------
    # For the captured and recorded gender individuals
    #--------------------------------------------------
    
    A1 = A2 = matrix(NA, numl, nG) # numl x nG
    for(ind in 1:numl)
    {
      if(sex.obs[ind] == 1)
      {
        U1[1:KnG] = log(pmatm.t[1:KnG]^rep(caphist2d[ind,], nG))
        U2[1:KnG] = log((1 - pmatm.t)[1:KnG] ^ rep(2*J-caphist2d[ind,], nG))
      }
      if(sex.obs[ind] == 0)
      {
        U1[1:KnG] = log(pmatf.t[1:KnG]^rep(caphist2d[ind,], nG))
        U2[1:KnG] = log((1 - pmatf.t)[1:KnG] ^ rep(2*J-caphist2d[ind,], nG))
      }
      A1[ind,] = colSums(U1, na.rm = T) # nG x 1 # Summing over the traps.
      A2[ind,] = colSums(U2, na.rm = T) # nG x 1 # Summing over the traps.
    }
    A1[A1 == -Inf] = -bignum; A2[A2 == -Inf] = -bignum
    part1 =  log(psi) + log(theta^sex.obs) + log((1-theta)^(1-sex.obs)) +
      apply(A1 + A2, 1, logMeanExp, na.rm = T) # numl x 1
    
    #--------------------------------------------------
    # For the captured but missing gender individuals
    #--------------------------------------------------
    if(numcap > numl){
      III = c(1:M)[yi00 > 0 & c(1:M) > numl] # (numcap - numl) x 1
      B1 = B2 = B3 = B4 = matrix(NA, (numcap - numl), nG) # (numcap - numl) x nG
      i = 0
      for(ind in III)
      {
        i = i + 1
        U1[1:KnG] = log(pmatm.t[1:KnG]^rep(caphist2d[ind,], nG))
        U2[1:KnG] = log((1 - pmatm.t)[1:KnG] ^ rep(2*J-caphist2d[ind,], nG))
        U3[1:KnG] = log(pmatf.t[1:KnG]^rep(caphist2d[ind,], nG))
        U4[1:KnG] = log((1 - pmatf.t)[1:KnG] ^ rep(2*J-caphist2d[ind,], nG))
        
        B1[i,] = colSums(U1, na.rm = T) # nG x 1 # Summing over the traps.
        B2[i,] = colSums(U2, na.rm = T) # nG x 1 # Summing over the traps.
        B3[i,] = colSums(U3, na.rm = T) # nG x 1 # Summing over the traps.
        B4[i,] = colSums(U4, na.rm = T) # nG x 1 # Summing over the traps.
      }
      
      B1[B1 == -Inf] = -bignum; B2[B2 == -Inf] = -bignum; 
      B3[B3 == -Inf] = -bignum; B4[B4 == -Inf] = -bignum; 
      
      part2 = log(psi) + log(theta * exp(apply(B1+B2, 1, logMeanExp, na.rm = T))+ 
                               (1-theta) * exp(apply(B3+B4, 1, logMeanExp, na.rm = T))) # (numcap - numl) x 1
      
    } # end of if(numcap > numl)
    
    #----------------------------------------------------
    # For the uncaptured and missing gender individuals
    #----------------------------------------------------
    U1[1:KnG] = log((1 - pmatm.t)[1:KnG] ^ rep(2*J, KnG))
    U2[1:KnG] = log((1 - pmatf.t)[1:KnG] ^ rep(2*J, KnG))
    
    C1 = colSums(U1, na.rm = T) # nG x 1 # Summing over the traps.
    C2 = colSums(U2, na.rm = T) # nG x 1 # Summing over the traps.
    C1[C1 == -Inf] = -bignum; C2[C2 == -Inf] = -bignum
    part3 = log(psi*(theta*mean(exp(C1),na.rm=T)+(1-theta)*mean(exp(C2),na.rm=T)) + (1-psi)) # scalar
    
    if(numcap > numl){loglik = sum(part1) + sum(part2) + (M-numcap) * part3}
    if(numcap == numl){loglik = sum(part1) + (M-numcap) * part3}
    return(loglik) 
    
  } # end of  if(Msexsigma == 1)
  #----------------------------------     
  if(Msexsigma == 0)
  {  
    logpmat.t = log(p0) - D.t * D.t / (2 * (sigma  ^ 2)) # ntrap x nG # Transpose of the usual pmat
    pmat.t = exp(logpmat.t)
    pmat.t[pmat.t < .Machine$double.xmin] = .Machine$double.xmin
    #-------------------------------
    # For the captured individuals
    #-------------------------------
    III = c(1:M)[yi00 > 0] # numcap x 1
    
    B1 = B2 = matrix(NA, numcap, nG) # numcap x nG
    U1 = U2 = matrix(NA, ntrap, nG) # ntrap x nG
    i = 0
    for(ind in III)
    {
      i = i+1
      U1[1:KnG] = log(pmat.t[1:KnG]^rep(caphist2d[ind,], nG))
      U2[1:KnG] = log((1 - pmat.t)[1:KnG] ^ rep(2*J-caphist2d[ind,], nG))
      
      B1[i,] = colSums(U1, na.rm = T) # nG x 1 # Summing over the traps.
      B2[i,] = colSums(U2, na.rm = T) # nG x 1 # Summing over the traps.
    }
    B1[B1 == -Inf] = -bignum; B2[B2 == -Inf] = -bignum;  
    part1 = log(psi) + apply(B1 + B2, 1, logMeanExp, na.rm = T) # numcap x 1
    
    #---------------------------------
    # For the uncaptured individuals
    #---------------------------------
    
    U1[1:KnG] = log((1 - pmat.t)[1:KnG] ^ rep(2*J, KnG))
    C1 = colSums(U1, na.rm = T) # nG x 1 # Summing over the traps.
    C1[C1 == -Inf] = -bignum; 
    part2 = log(psi * mean(exp(C1), na.rm = T) + (1-psi)) # scalar
    #----------------------------------     
    loglik = sum(part1) + (M-numcap) * part2
    return(loglik) 
    
  } # end of  if(Msexsigma == 0)
}

sim.partial.data = function(N, N.Male = NA, phi, omega0, sigma, Msexsigma, nrow_trap, ncol_trap, nocc, xlim, ylim, buffer)
{
  
  e2dist = function(x,y) # x is matrix mx2, y is matrix nx2
  {
    ind = sort(rep(1:nrow(y), nrow(x))) # [rep(1, nrow(x)), rep(2, nrow(x)), ..., rep(nrow(y), nrow(x))]'
    dvec = sqrt((x[,1] - y[ind,1]) ^ 2 + (x[,2] - y[ind,2]) ^ 2) # [d(1,1), d(2,1), ...,d(nrow(x), 1), ..., d(1,nrow(y)), d(2,nrow(y)),...,d(nrow(x), nrow(y))]
    dmat = matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F) 
    return(dmat)
  }
  
  # simulate a population of activity centers 
  act.center = cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
  act.x = act.center[,1]
  act.y = act.center[,2]
  state.space.area = (xlim[2] - xlim[1])*(ylim[2] - ylim[1])
  trap.locations = as.matrix(expand.grid(seq(xlim[1] + buffer, xlim[2] - buffer, length.out = ncol_trap), 
                                         seq(ylim[1] + buffer, ylim[2] - buffer, length.out = nrow_trap)))
  traploc.x = trap.locations[,1]
  traploc.y = trap.locations[,2]
  ntrap = nrow_trap * ncol_trap
  
  if(Msexsigma == 1)
  {
    # Getting the genders of each individual
    ind = sample(1:N, N.Male, replace = F)
    sex.true=rep(0,N)
    sex.true[ind]=1
    sig = rep(sigma[1], N)
    sig[sex.true == 0] = sigma[2]
  }
  
  if(Msexsigma == 0){ sig = sigma }
  
  # Trap entry probability
  logpimat = log(omega0) - e2dist(act.center, trap.locations) ^ 2 / (2 * (sig ^ 2)) # Division is always column wise 
  
  # Getting the captures
  left = right = array(0,dim = c(N,ntrap,nocc))
  for (i in 1:N){ for (j in 1:ntrap){ for (k in 1:nocc){
    u = runif(1,0,1)
    if (log(u) <= logpimat[i,j])
    {
      left[i,j,k] = rbinom(1,1,phi)
      right[i,j,k] = rbinom(1,1,phi)
    } 
  }}}
  
  # Put known IDfixed guys first.
  known.ID = apply(left+right, 1, function(amat){sum(amat == 2) > 0}) # individuals with at least one simultaneous capture (l,r)=(1,1)
  known.ID.indx = c(1:N)[known.ID] # IDs of fully identified individuals
  IDfixed = sum(known.ID) # No. of fully identified individuals
  
  kp.left = apply(left,1,sum) > 0   # Indicator of individuals captured on left side
  kp.right = apply(right,1,sum) > 0 # Indicator of individuals captured on right side
  
  numl = sum(kp.left)
  numr = sum(kp.right)
  if(numl == 0 & numr > 0){warning("Zero number of captured individuals by detector 1")}
  if(numl > 0 & numr == 0){warning("Zero number of captured individuals by detector 2")}
  if(numl == 0 & numr == 0){warning("Zero number of captured individuals by both detector 1 and detector 2")}
  
  if (IDfixed == 0) known = 'none'
  if (IDfixed > 0)
  {
    known = 'some'
    logic1 = all.equal(kp.left, known.ID)
    logic2 = all.equal(kp.right, known.ID)
    if (logic1 == T & logic2 == T) known = 'ALL'
  }
  
  original.left = left
  original.right = right
  
  if (known != 'ALL')
  {
    which.left.obsd = c(1:N)[kp.left] 
    left = left[c(known.ID.indx, setdiff(which.left.obsd, known.ID.indx)),,]
    if(Msexsigma == 1){ sexl = sex.true[c(known.ID.indx, setdiff(which.left.obsd, known.ID.indx))] }
    if(Msexsigma == 0){ sexl = NULL }
    which.left.obsd = c(known.ID.indx, setdiff(which.left.obsd, known.ID.indx))
    
    which.right.obsd = c(1:N)[kp.right] 
    right = right[c(known.ID.indx, setdiff(which.right.obsd, known.ID.indx)),,]
    if(Msexsigma == 1){ sexr = sex.true[c(known.ID.indx, setdiff(which.right.obsd, known.ID.indx))] }
    if(Msexsigma == 0){ sexr = NULL }
    which.right.obsd = c(known.ID.indx, setdiff(which.right.obsd, known.ID.indx))
  }
  
  if (known == 'ALL')
  {
    cap = apply(left+right,1,sum)>0 
    left = left[cap,,]
    right = right[cap,,]
    if(Msexsigma == 1){ sexl = sexr = sex.true[cap] }
    if(Msexsigma == 0)
    {
      N.Male = NA
      sexl = sexr = NULL
      sex.true = NULL
      logsigmaf = NA
    }
    numl = numr = sum(cap)
    which.left.obsd = which.right.obsd = NULL
  }
  
  # Trap deployment file with nocc + 3 columns - TrapID, Trap location (x), Trap location (y), Trap activity indicator for each occasion
  tdf = cbind(1:ntrap, trap.locations, matrix(1, ntrap, nocc)); mask = tdf[, 4:(nocc+3)]
  dimnames(tdf)[[2]] = c('TrapID', 'Easting', 'Northing', 1:nocc)
  
  left_edf = right_edf = NULL
  if(numl > 0){
    if(numl == 1){left = array(left, c(1, ntrap, nocc))}
    for(i in 1:numl){ for(k in 1:nocc){ for(j in 1:ntrap){ 
      if(left[i,j,k] == 1 & mask[j,k] == 1) left_edf = rbind(left_edf, c(1, i, k, j))
    }}}
    diml1 = dim(left_edf)[1]
    left_edf = cbind(left_edf, rep(12, diml1))
    if(IDfixed < numl){
      ind1 = c(1:diml1)[left_edf[,2]>IDfixed]
      left_edf[ind1, 5] = 1
    } 
    dimnames(left_edf)[2] = list(c('Session', 'Individual', 'Occasion', 'Trap', 'Status')) # `Status' means which detector this data corresponds to
  }
  if(numr > 0){
    if(numr == 1){right = array(right, c(1, ntrap, nocc))}
    for(i in 1:numr){ for(k in 1:nocc){ for(j in 1:ntrap){ 
      if(right[i,j,k] == 1 & mask[j,k] == 1) right_edf = rbind(right_edf, c(1, i, k, j))
    }}} 
    dimr1 = dim(right_edf)[1]
    right_edf = cbind(right_edf, rep(12, dimr1))
    if(IDfixed < numr){
      ind2 = c(1:dimr1)[right_edf[,2]>IDfixed]
      right_edf[ind2, 5] = 2
    } 
  dimnames(right_edf)[2] = list(c('Session', 'Individual', 'Occasion', 'Trap', 'Status')) # `Status' means which detector this data corresponds to
  }
  
  # Encounter deployment files with 5 columns : Session, Individual, Occasion, TrapID, Status
  
  if(Msexsigma == 1){
    if(numl > 0 & numr > 0){
      if(numl > numr){ sex.mat = cbind(sexl[1:numl], c(sexr[1:numr], rep(NA, numl - numr)))}  
      if(numl < numr){ sex.mat = cbind(c(sexl[1:numl], rep(NA, numr - numl)), sexr[1:numr])}
      if(numl == numr){ sex.mat = cbind(sexl[1:numl], sexr[1:numr])}
      dimnames(sex.mat)[[2]] = c('sex1', 'sex2')
      sex.data = matrix(NA, max(numl, numr), 4) 
      sex.data[1:numl, 1] = 1:numl; sex.data[1:numr, 3] = 1:numr
      sex.data[1:numl, 2] = 'Male'; sex.data[1:numr, 4] = 'Male'
      sex.data[which(sex.mat[,1] == 0, arr.ind = T), 2] = 'Female'
      sex.data[which(sex.mat[,2] == 0, arr.ind = T), 4] = 'Female'
      dimnames(sex.data)[[2]] = c('Individual1', 'sex1', 'Individual2', 'sex2')
    }
    if(numl > 0 & numr == 0){
      sex.mat = cbind(sexl[1:numl], rep(0, numl))   
      dimnames(sex.mat)[[2]] = c('sex1', 'sex2')
      sex.data = matrix(NA, numl, 4) 
      sex.data[1:numl, 1] = 1:numl
      sex.data[1:numl, 2] = 'Male'
      sex.data[which(sex.mat[,1] == 0, arr.ind = T), 2] = 'Female'
      dimnames(sex.data)[[2]] = c('Individual1', 'sex1', 'Individual2', 'sex2')
    }
    if(numl == 0 & numr > 0){
      sex.mat = cbind(rep(0, numr), sexr[1:numr])
      dimnames(sex.mat)[[2]] = c('sex1', 'sex2')
      sex.data = matrix(NA, numr, 4) 
      sex.data[1:numr, 3] = 1:numr
      sex.data[1:numr, 4] = 'Male'
      sex.data[which(sex.mat[,2] == 0, arr.ind = T), 4] = 'Female'
      dimnames(sex.data)[[2]] = c('Individual1', 'sex1', 'Individual2', 'sex2')
    }
    if(numl == 0 & numr == 0){ sex.mat = sex.data = NULL }
  }
  
  if(Msexsigma == 0){  sex.mat = sex.data = NULL }
  
  out = list(edf1 = left_edf, # matrix
             edf2 = right_edf, # matrix
             sex.data = sex.data, # matrix
             tdf = tdf, # matrix
             xlim = xlim, #/ scale # 2 x 1
             ylim = ylim, #/ scale # 2 x 1
             buffer = buffer, #/ scale # scalar
             Msexsigma.given = Msexsigma, # binary
             N.given = N, # scalar
             N.Male.given = N.Male, # scalar
             phi.given = phi, # scalar
             omega0.given = omega0, # scalar
             sigma.given = sigma # vector 2x1 if Msexsigma = 1, scalar if Msexsigma = 0
  )
  return(out)
  
} # end of sim.partial.data

sim.partial.data.ppl = function(phi, basepar, sigma, sex, z, L, S, trap.locations, nocc, model)
{
  e2dist = function(x,y) # x is matrix mx2, y is matrix nx2
  {
    ind = sort(rep(1:nrow(y), nrow(x))) # [rep(1, nrow(x)), rep(2, nrow(x)), ..., rep(nrow(y), nrow(x))]'
    dvec = sqrt((x[,1] - y[ind,1]) ^ 2 + (x[,2] - y[ind,2]) ^ 2) # [d(1,1), d(2,1), ...,d(nrow(x), 1), ..., d(1,nrow(y)), d(2,nrow(y)),...,d(nrow(x), nrow(y))]
    dmat = matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F) 
    return(dmat)
  }
  
  M = dim(S)[1]
  ntrap = dim(trap.locations)[1] # nrow_trap * ncol_trap
  if(model %in% c(1,2)){ sig = rep(sigma[1], M); sig[sex == 0] = sigma[2]} # Msexsigma = 1
  if(model %in% c(3,4)){ sig = sigma; sexl = NULL; sexl.cap = NULL } # Msexsigma = 0
  logpimat = log(basepar) - e2dist(S, trap.locations) ^ 2 / (2 * (sig ^ 2)) # Division is always column wise 
  left = right.star = array(0, dim = c(M,ntrap,nocc))
  #==========================================
  # Dey et al. (2019)
  if(model %in% c(1,3)){
    for (i in 1:M){ for (j in 1:ntrap){ for (k in 1:nocc){
      u = runif(1,0,1)
      if (log(u) <= logpimat[i,j]){
          left[i,j,k] = rbinom(1, 1, phi); right.star[i,j,k] = rbinom(1, 1, phi)
      }  }}} }
  #==========================================
  # Royle (2015)
  if(model %in% c(2,4)){
    for (k in 1:nocc){ 
     left[,,k] = matrix(rbinom(M*ntrap, 1, exp(logpimat)), M, ntrap);
     right.star[,,k] = matrix(rbinom(M*ntrap, 1, exp(logpimat)), M, ntrap)
    }  }
  #==========================================
  left[z==0,,] = 0
  right.star[z==0,,] = 0
  right = right.star[L,,]
  out = list(crdata1 = left, crdata2 = right)
  return(out)
  
} # end of sim.gender.data.partial


#====================================================================
# compute Monte Carlo standard error of a functional (FUNC) of a Markov chain (x) using the Subsampling Bootstrap Method (SBM)
mcse = function(x, FUNC,  batchsize = floor(sqrt(length(x)))) {
  n = length(x)
  nb = n - batchsize + 1
  bval = rep(NA,nb)
  for (j in 1:nb) {
    ind = j:(batchsize + j - 1)
    bval[j] = FUNC(x[ind])
  }
  var.bval = var(bval)*(nb - 1) * batchsize / nb
  return(list(se = sqrt(var.bval / n), batchsize = batchsize))
}

#====================================================================
# Good initial choice of `L' vector.
Lvec.fn =  function(left, right, numl, numr, sexl, sexr, traploc, IDfixed, nloopL, Msexsigma)
{
  e2dist = function(x,y)
  {
    ind = sort(rep(1:nrow(y), nrow(x))) # [rep(1, nrow(x)), rep(2, nrow(x)), ..., rep(nrow(y), nrow(x))]'
    dvec = sqrt((x[,1] - y[ind,1]) ^ 2 + (x[,2] - y[ind,2]) ^ 2) # [d(1,1), d(2,1), ...,d(nrow(x), 1), ..., d(1,nrow(y)), d(2,nrow(y)),...,d(nrow(x), nrow(y))]
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F) 
  }

  # function needs to spit out initial L and activity centers
  M = dim(left)[1]
  ntrap = dim(left)[2]
  
  idkp = (IDfixed + 1):M
 
  numl = numl - IDfixed
  numr = numr - IDfixed
  
  left2d = apply(left[idkp,,], c(1,2), sum, na.rm = T) # total left captures per trap for indinvidals (IDfixed+1):M
  right2d = apply(right[idkp,,],c(1,2), sum, na.rm = T) # total right captures per trap for indinvidals (IDfixed+1):M
  if(Msexsigma == 1)
  {
    sexl = sexl[idkp]
    sexr = sexr[idkp]
    
    temp.lm = which(sexl == 1, arr.ind = T)
    temp.lf = which(sexl == 0, arr.ind = T)
    
    temp.rm = which(sexr == 1, arr.ind = T)
    temp.rf = which(sexr == 0, arr.ind = T)
    
    temp.lr = which(sexl == sexr, arr.ind = T)
  }
  traploc = as.matrix(traploc)
  sbar.left = matrix(NA,nrow = (M - IDfixed),ncol = 2)
  sbar.right = matrix(NA,nrow = (M - IDfixed),ncol = 2)
  
  #------------------------------------------------------------------------
  # Getting an initial value for activity centres sbar.left and sbar.right
  #------------------------------------------------------------------------
  for (i in 1:(M - IDfixed))
  {
    #-----------  
    # FOR LSIs
    #-----------
    if (sum(left2d[i,]) > 0) # should always be satisfied i.e LSi i got captured at least once or not
    {  
      # get the traps where ith LSI got captured
      traps.loc = matrix(traploc[left2d[i,] > 0,], ncol = 2, byrow = F)
      sbar.left[i,] = apply(traps.loc, 2, weighted.mean, w = left2d[i,][left2d[i,] > 0]) #c(mean(trps[,1]),mean(trps[,2]))
    }
    if (sum(left2d[i,]) == 0) sbar.left[i,] = traploc[sample(1:ntrap,1),] # if LSi i never got captured
    
    #-----------  
    # FOR RSIs
    #-----------
    
    if (sum(right2d[i,]) > 0) # should always be satisfied i.e RSi i got captured at least once or not
    { 
      # get the traps where ith LSI got captured
      traps.loc = matrix(traploc[right2d[i,] > 0,], ncol = 2, byrow = F)
      sbar.right[i,] = apply(traps.loc, 2, weighted.mean, w = right2d[i,][right2d[i,] > 0]) #c(mean(trps[,1]),mean(trps[,2]))
    }
    if (sum(right2d[i,]) == 0) sbar.right[i,] = traploc[sample(1:ntrap,1),] # if RSi i never got captured
    
  } # end of for(i in 1:(M-IDfixed))
  
  #Getting euclidean distance between each right act centre sbar.right and each leftt act centre sbar.left
  D = e2dist(sbar.right,sbar.left) # dimension (M-IDfixed) x (M-IDfixed)
  
  #---------------------------------------------------------
  # Initializing L i.e linking RSI with LSI wrt same sex
  #---------------------------------------------------------

  L = sample(1:(M - IDfixed), (M - IDfixed), replace = F)
  
  Q = sum( D[cbind(1:(M - IDfixed),L)]  ) # sum of all the distances between pair sbar.right & sbar.left
  
  for (loop in 1:nloopL)
  {
    for (i in 1:(M - IDfixed))
    {
      # Checking distance sum by linking RSI i to other LSI other than L[i]
      curr.spot = L[i]
      Qtmp = c()  
      if(Msexsigma == 1)
      {
        if(sexr[i] == 1){ ind1 = setdiff(temp.lm, i) }
        if(sexr[i] == 0){ ind1 = setdiff(temp.lf, i) }
      }
      if(Msexsigma == 0){  ind1 = setdiff(1:(M - IDfixed), i) }
     
      for (k in ind1)
      {
        L[i] = L[k];  L[k] = curr.spot # swap
        Qtmp = c(Qtmp, sum( D[cbind(1:(M - IDfixed),L)] ))
        L[k] = L[i];  L[i] = curr.spot # swap back
      } # end of for(k in ind1)
      
      if (min(Qtmp, na.rm = T) < Q )
      {
        # Make the swap
        which_L = ind1[Qtmp == min(Qtmp)][1] # ind1[which(Qtmp == min(Qtmp), arr.ind = T)][1]
        L[i] = L[which_L]
        L[which_L] = curr.spot
        
        Q = min(Qtmp, na.rm = T)
      } # end of if(min(Qtmp, na.rm = T) < Q )
   
    } # close loop over 1:(M-IDfixed) [ =length of L ]
  } # close loop over 1:50
  
  if (IDfixed > 0) L = c((1:IDfixed) ,  IDfixed + L )
  if (IDfixed == 0) L = L
  
  return(L)
}

SCR23darray = function (edf, tdf)  # same as 'SCR23darray' function SCRbayes package
{
  uq <- paste(edf[, 2], edf[, 3], edf[, 4])
  uq <- unique(uq)
  if (any(table(uq) > 1)) 
    cat("Duplicate captures (same individual, trap, occasion) present in data set, these are not used", 
        fill = TRUE)
  nind <- max(edf[, 2])
  ntraps <- nrow(tdf)
  nperiods <- ncol(tdf) - 3
  per.id <- as.numeric(dimnames(tdf)[[2]][4:ncol(tdf)])
  ind.id <- edf[, 2]
  trap.id <- edf[, 4]
  if (length(per.id) != length(min(per.id):max(per.id))) {
    x <- 1:nperiods
    names(x) <- as.character(per.id)
    per.id <- x[as.character(edf[, 3])]
  }
  else {
    per.id <- edf[, 3]
  }
  y <- array(0, c(nind, ntraps, nperiods))
  tmp <- cbind(ind.id, trap.id, per.id)
  y[tmp] <- 1
  y
}