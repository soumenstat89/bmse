# rm(list = ls())
options(digits = 8)

start.time = Sys.time()
ts = format(Sys.time(), "%d%m%y_%H%M%S")


fff0 = list.files(path = getwd(), pattern = "Bilateral_MCMC")
folderNames1 = c()

for(model in 1:4) { folderNames1 = c(folderNames1, paste(getwd(), '/', fff0[grep(paste("_M", model, "_", sep = ""), fff0)], '', sep = ''))}

for(model in 1:4){
  cat('Model = ', model, '\n', sep = '')
  folderpath = folderNames1[model]
  
  #----------------------------
  # Handling of simulated data
  #----------------------------
  
  info1 = read.csv(paste(folderpath, '/info1.csv', sep = ''), sep = ',', header = T)
  Msexsigma = info1[, 'Msexsigma'] # set at MCMC sampling
  xlim = unlist(info1[, c('xlim.1', 'xlim.2')])
  ylim = unlist(info1[, c('ylim.1', 'ylim.2')])
  buffer = info1[, 'buffer']
  ndraws = info1[, 'ndraws']
  burnin = info1[, 'burnin']
  K = info1[, 'ntrap']
  J = info1[, 'nocc']
  M = info1[, 'M']
  
  numl = info1[,'num1']
  numr = info1[,'num2']
  IDfixed = info1[, 'IDfixed']
  known = info1[, 'known']
  
  if(Msexsigma == 1)
  {
    sex.data = read.csv(paste(folderpath, '/sex_captured.csv', sep = ''), sep = ',', header = T) # max(numl, numr) x 1
    sexl.data = na.omit(sex.data[,'sex1']) # numl x 1
    sexr.data = na.omit(sex.data[,'sex2']) # numr x 1
    sexl.obs = rep(1, numl)
    sexr.obs = rep(1, numr)
    if(sum(sexl.data == 'Female') > 0){ sexl.obs[sexl.data == 'Female'] = 0} # numl x 1
    if(sum(sexr.data == 'Female') > 0){ sexr.obs[sexr.data == 'Female'] = 0} # numr x 1
  } 
  if(Msexsigma == 0)
  {  
    sexl.obs = sexr.obs = NULL
  }
  
  
  #------------------------------------------------------  
  #--------------------------
  # Trap locations
  #--------------------------
  tdf = read.csv(paste(folderpath, '/trap_deployment_file.csv', sep = ''), sep = ',', header = T)
  dimnames(tdf)[[2]] = c('TrapID', 'Easting', 'Northing', 1:J)
  traploc = as.matrix(tdf[,c('Easting', 'Northing')] )
  mask = tdf[, 4:(J+3)]
  
  
  #------------------------------------------------------  
  #--------------------------
  # Capture-recapture data
  #--------------------------
  left_edf = read.csv(paste(folderpath, '/encounter_data_file1.csv', sep = ''), sep = ',', header = T)
  right_edf = read.csv(paste(folderpath, '/encounter_data_file2.csv', sep = ''), sep = ',', header = T)
  
  if(   length(unique(left_edf[,2])) != length(min(left_edf[,2]):max(left_edf[,2])) ) {
    cat("Error: individuals not numbered sequentially, renumbering them now",fill=TRUE)
    left_edf[,2]<- as.numeric(factor(left_edf[,2]))
  }
  
  if(   length(unique(right_edf[,2])) != length(min(right_edf[,2]):max(right_edf[,2])) ) {
    cat("Error: individuals not numbered sequentially, renumbering them now",fill=TRUE)
    right_edf[,2]<- as.numeric(factor(right_edf[,2]))
  }
  
  left = left.obs = SCR23darray(left_edf[,c('Session', 'Individual', 'Occasion', 'Trap')], tdf) # numl x K x J # left3d = left
  right = right.obs = SCR23darray(right_edf[,c('Session', 'Individual', 'Occasion', 'Trap')], tdf) # numr x K x J # right3d = right
  
  # Also re-label data sets if needed so that the left side always has more individuals
  if (known!="ALL" & numl < numr)
  { 
    a = left; left = right; right = a
    b = numl; numl = numr; numr = b
    if(Msexsigma == 1)
    {
      sexl = sexr.obs # updated numl x 1
      sexr = sexl.obs # updated numr x 1
    }   
  }
  
  # Using the mask to make some false-positives zero
  for(t in 1:J){ for(k in 1:K){ 
    if(mask[k, t] == 0) {left[1:numl, k, t] = 0; right[1:numr, k, t] = 0}
  }}
  
  left = abind(left, array(0, dim = c( M - numl, K, J)), along = 1)
  right = abind(right, array(0, dim = c( M - numr, K, J)), along = 1)
  
  left2d = apply(left, c(1,2), sum) # no. of times individual i got captured at trap k on LEFT side
  right2d = apply(right, c(1,2), sum) # no. of times individual i got captured at trap k on RIGHT side
  
  yli00 = rowSums(left); 
  yri00 = rowSums(right);
  #----------------------------------------------
  
  post = read.csv(paste0(folderpath, '/markovchain.txt', sep = ''), sep = ',', header = T)
  col.nms = dimnames(post)[[2]]
  tot.length = (ndraws - burnin)
   
  post = matrix(unlist(post[(burnin + 1):ndraws,]), tot.length, dim(post)[2])
  dimnames(post)[[2]] = col.nms
  
  post.N = post[, 'N']; 
  post.psi = post[, 'psi']; post.logitpsi = logit(post.psi)
  
  if(model == 1 | model == 3)
  {
    post.phi = post[, 'phi']; post.logitphi = logit(post.phi)
    post.omega0 = post[, 'omega0']; post.logitomega0 = logit(post.omega0)
  }
  if(model == 2 | model == 4)
  {
    post.p0 = post[, 'p0']; post.logitp0 = logit(post.p0)
  }
  if(model == 1 | model == 2)
  { 
    post.N.Male = post[, 'N.Male']; 
    post.theta = post[, 'theta']; post.logittheta = logit(post.theta)
    post.sigmam = post[, 'sigmam']; post.logsigmam = log(post.sigmam)
    post.sigmaf = post[, 'sigmaf']; post.logsigmaf = log(post.sigmaf)
    post.sex = read.csv(paste0(folderpath, '/markovchain.sex.txt', sep = ''), sep = ',', header = T)[(burnin + 1):ndraws,]
    post.sex.miss = post.sex[, (numl+1):M]
    sex.obs = post.sex[1,1:numl] 
    
  }
  if(model == 3 | model == 4)
  {
    post.sigma = post[, 'sigma']; post.logsigma = log(post.sigma)
  }
  
  post.sx = read.csv(paste0(folderpath, '/markovchain.sx.txt', sep = ''), sep = ',', header = T)[(burnin + 1):ndraws,]
  post.sx = matrix(unlist(post.sx), tot.length, M)
  post.sy = read.csv(paste0(folderpath, '/markovchain.sy.txt', sep = ''), sep = ',', header = T)[(burnin + 1):ndraws,]
  post.sy = matrix(unlist(post.sy), tot.length, M)
  post.s = cbind(post.sx, post.sy)
  
  if (known != 'ALL'){
    post.L = read.csv(paste0(folderpath, '/markovchain.L.txt', sep = ''), sep = ',', header = T)[(burnin + 1):ndraws,]
    post.L = matrix(unlist(post.L), tot.length, M)
  }
  
  if (known == 'ALL'){  post.L = matrix(1:M, nrow = ndraws - burnin, ncol = M, byrow = T)}
  
  post.z = read.csv(paste0(folderpath, '/markovchain.z.txt', sep = ''), sep = ',', header = T)[(burnin + 1):ndraws,]
  post.z = matrix(unlist(post.z), tot.length, M)
  
#==========
# Data
#==========
R = sqrt((xlim[2]-xlim[1])^2 + (ylim[2]-ylim[1])^2) # R = Maximum stretch of the state space

#======================
# Log-likelihood chain
#======================

loglik.u0zs.chain = unlist(read.csv(paste(folderpath, '/markovchain.loglikelihood.u0zs.txt', sep = ''), sep = ',', header = T)) # tot.length x 1

if(model == 1 | model  == 2)
{
  logfactor.sex = log(post.theta^(post.z*post.sex)) + log((1 - post.theta)^(1 - post.z*post.sex)) # tot.length x M
}
#======================
# Log-prior chain
#======================
# Prior (Transformed) 
# <[Priors are transformed because tuning density uses transformed parameters]>
logprior.logitpsi = log((exp(post.logitpsi))/((1 + exp(post.logitpsi))^2))
logprior.L =  - lgamma(M+1) # = - log(factorial(M)) = log(1/M!) 
if(model == 1 | model == 3){
  logprior.logitphi = log((exp(post.logitphi))/((1 + exp(post.logitphi))^2)) 
  logprior.logitomega0 = log((exp(post.logitomega0))/((1 + exp(post.logitomega0))^2))
}
if(model == 2 | model == 4){
  logprior.logitp0 = log((exp(post.logitp0))/((1 + exp(post.logitp0))^2))
}
if(model == 1 | model == 2)
{ 
  logprior.logittheta = log((exp(post.logittheta))/((1 + exp(post.logittheta))^2))
  logprior.logsigmam = log((exp(post.logsigmam))/R) 
  logprior.logsigmaf = log((exp(post.logsigmaf))/R)
}
if(model == 3 | model == 4){ logprior.logsigma = log((exp(post.logsigma))/R) }

if(model == 1){ logprior.chain = logprior.logitpsi + logprior.logittheta + logprior.logitphi + logprior.logitomega0 + 
                                   logprior.logsigmam + logprior.logsigmaf + logprior.L  }
if(model == 2){ logprior.chain = logprior.logitpsi + logprior.logittheta + logprior.logitp0 + 
                                   logprior.logsigmam + logprior.logsigmaf + logprior.L }
if(model == 3){ logprior.chain = logprior.logitpsi + logprior.logitphi + logprior.logitomega0 + 
                                   logprior.logsigma + logprior.L }
if(model == 4){ logprior.chain = logprior.logitpsi + logprior.logitp0 + 
                                   logprior.logsigma + logprior.L }

if(model == 1){ post.chains = cbind(post.logitpsi, post.logittheta, post.logitphi, post.logitomega0, post.logsigmam, post.logsigmaf)}
if(model == 2){ post.chains = cbind(post.logitpsi, post.logittheta, post.logitp0, post.logsigmam, post.logsigmaf) }
if(model == 3){ post.chains = cbind(post.logitpsi, post.logitphi, post.logitomega0, post.logsigma) }
if(model == 4){ post.chains = cbind(post.logitpsi, post.logitp0, post.logsigma) }


#================================================
# Parameters for tuning function Multivariate-t (CHOICE - 2)
#================================================

mu.est = apply(post.chains, 2, mean, na.rm = T)
Sigma.est.choice2 = ((dim(post.chains)[1]-1)/dim(post.chains)[1])*cov(post.chains)
diag(Sigma.est.choice2)[diag(Sigma.est.choice2)< 10^(-5)] = 10^(-5)
out2 = log.g.chain2 = NULL


#==============================================
# Gelfand-Dey marginal likelihood estimator
#==============================================
# g(L) = pi(L) i.e Takng prior of L as the tuning density of L 
log.g.norm.chain = dmvnorm(x = post.chains, mean = mu.est,
                            sigma = Sigma.est.choice2, log = T) + 
                    logprior.L
log.g.chain2 = cbind(log.g.chain2, log.g.norm.chain)
logh.chain = log.g.norm.chain-(loglik.u0zs.chain+logprior.chain)
logmarglik.normalg = gdmean(logh.chain)
out2 = rbind(out2, c(logmarglik.normalg,NA))

degf = c(10, 100, 500, 1000, 10000)

for(i in 1:length(degf))
{
  
  log.g.t1.chain = dmvt(x = post.chains, delta = mu.est, sigma = Sigma.est.choice2,
                        df =degf[i], log = T, type = "shifted") +
                    logprior.L  
  log.g.chain2 = cbind(log.g.chain2, log.g.t1.chain)
  logh.chain = log.g.t1.chain-(loglik.u0zs.chain+logprior.chain)
  logmarglik.tg1 = gdmean(logh.chain)
  out2 = rbind(out2, c(logmarglik.tg1,NA))
  
}

conflev = c(0.9, 0.95, 0.99)
q.chisq = qchisq(conflev, df = dim(post.chains)[2], ncp = 0, lower.tail = T)
scoreval = apply(post.chains, 1, scorefn, mu = mu.est, Sigma = Sigma.est.choice2)
for(lev in 1:length(conflev))
{
  indices = scoreval < q.chisq[lev]
  log.g.trunc.norm.chain = dmvtruncnorm(x = post.chains[indices,], mu = mu.est,
                                        Sigma = Sigma.est.choice2,
                                        conf.lev = conflev[lev], logg = T) + 
                           logprior.L
  if(length(log.g.trunc.norm.chain) == tot.length) {log.g.chain2 = cbind(log.g.chain2, log.g.trunc.norm.chain)}
  if(length(log.g.trunc.norm.chain) < tot.length) {log.g.chain2 = cbind(log.g.chain2, c(log.g.trunc.norm.chain, rep(NA, tot.length - length(log.g.trunc.norm.chain))))}
  logh.chain = log.g.trunc.norm.chain-(loglik.u0zs.chain[indices]+logprior.chain[indices])
  logmarglik.gewekeg = gdmean(logh.chain)
  out2 = rbind(out2, c(logmarglik.gewekeg, NA))
  
}

tuning_den_names = c('GD.Normal',  paste('GD.t', c(10, 100, 500, 1000, 10000), sep = ''),
                     paste('GD.Trunc_normal_conf', c(0.9,0.95, 0.99), sep = ''))

dimnames(out2)= list(c(tuning_den_names), c('COL1', 'COL2'))
dimnames(log.g.chain2)[[2]] = tuning_den_names
write.csv(out2, file = paste(folderpath, '/marginal.IL.csv', sep = ''), quote = F, row.names = T)
write.csv(log.g.chain2, file = paste(folderpath, '/chain.log.g.IL.txt', sep = ''), quote = F,row.names = F)
write.csv(logprior.chain, file = paste(folderpath, '/logLik_logPrior.chain.IL.txt', sep = ''), quote = F,row.names = F)


} # end of for(model in 1:4)

#====================================================
save.image(paste0(folderpath, '/savingRimage_GDBF_IL.RData', sep = ''))
time.taken = Sys.time() - start.time; print(time.taken)
#=================================================
