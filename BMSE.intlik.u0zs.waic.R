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
  state.space.area = (xlim[2] - xlim[1])*(ylim[2] - ylim[1]) 
  # Getting the possible activity centers after discretising the state space (0, 5) x (0, 7)
  sgrid.dx = sgrid.dy = 0.1 
  R = sqrt((xlim[2]-xlim[1])^2 + (ylim[2]-ylim[1])^2) # R = Maximum stretch of the state space
  sgrid.nx = round((xlim[2]-xlim[1])/sgrid.dx); sgrid.ny = round((ylim[2]-ylim[1])/sgrid.dy)
  sgridloc = as.matrix(expand.grid(seq(xlim[1]+sgrid.dx/2, xlim[2]-sgrid.dx/2, by = sgrid.dx),
                                   seq(ylim[1]+sgrid.dy/2, ylim[2]-sgrid.dy/2, by = sgrid.dy)))
  nG = dim(sgridloc)[1] 
  if(model == 1 | model == 2)
  {
    data1 = list(left = left, left2d = left2d, right = right, right2d = right2d, 
                 numl = numl, sex.obs = sex.obs, traploc = traploc,  
                 J = J, M = M, sgridloc = sgridloc, state.space.area = state.space.area)
  }
  
  if(model == 3 | model == 4)
  {
    data1 = list(left = left, left2d = left2d, right = right, right2d = right2d, 
                 numl = numl, sex.obs = NULL, traploc = traploc,  
                 J = J, M = M, sgridloc = sgridloc, state.space.area = state.space.area)
  }
  
#============================================================
# Function to compute log-marginal likelihood over x0, z, s
#============================================================
 
if(model == 1){
simi.u0zs = function(draw)
{
  library(matrixStats)
  chain.val = list(psi = post.psi[draw], theta = post.theta[draw], 
                 phi = post.phi[draw], omega0 = post.omega0[draw], 
                 sigmam = post.sigmam[draw], sigmaf = post.sigmaf[draw], 
                 L = post.L[draw,], Msexsigma = Msexsigma)
  loglik = logLfn.piscr.u0zs (chain.val, data1)
  return(loglik)
}
}

if(model == 2){
simi.u0zs = function(draw)
{
  library(matrixStats)
  chain.val = list(psi = post.psi[draw], theta = post.theta[draw], 
                   p0 = post.p0[draw], 
                   sigmam = post.sigmam[draw], sigmaf = post.sigmaf[draw], 
                   L = post.L[draw,], Msexsigma = Msexsigma)
  loglik = logLfn.rpiscr.u0zs (chain.val, data1)
  return(loglik)
}
}


if(model == 3){
  simi.u0zs = function(draw)
  {
    library(matrixStats)
    chain.val = list(psi = post.psi[draw], 
                     phi = post.phi[draw], omega0 = post.omega0[draw], 
                     sigma = post.sigma[draw], 
                     L = post.L[draw,], Msexsigma = Msexsigma)
    loglik = logLfn.piscr.u0zs (chain.val, data1)
    return(loglik)
  }
}

if(model == 4){
  simi.u0zs = function(draw)
  {
    library(matrixStats)
    chain.val = list(psi = post.psi[draw], 
                     p0 = post.p0[draw], 
                     sigma = post.sigma[draw],
                     L = post.L[draw,], Msexsigma = Msexsigma)
    loglik = logLfn.rpiscr.u0zs (chain.val, data1)
    return(loglik)
  }
}
#=============================================================
  #===============================================================================
  # Function to compute log-likelihood at each MCMC iteration of every data point
  #===============================================================================
  
  if(model == 1 | model == 3){
    
    simi.waic = function(draw){
      library(matrixStats)
      right.star = right[order(post.L[draw,]),,]
      right2d.star = right2d[order(post.L[draw,]),]
      n0mat = apply(left + right.star, c(1,2), function(a){sum(a > 0)}) # trapwise count of recorded entries for each of the M individual
      n0vec = rowSums(n0mat) # giving total count of recorded entries for each of the M individual
      yli00 = rowSums(left2d); yri00.star = rowSums(right2d.star); yi00 = yli00 + yri00.star 
      if(model == 1)
      {
        
        logpimat = logpimatfn(logitbasepar = post.logitomega0[draw], s.left = cbind(c(post.sx[draw,]), c(post.sy[draw,])),
                        traploc = traploc, sex = post.sex[draw,], Msexsigma = Msexsigma,
                        logsigmam = post.logsigmam[draw], logsigmaf = post.logsigmaf[draw])
      }
      if(model == 3){
        
        logpimat = logpimatfn(logitbasepar = post.logitomega0[draw], s.left = cbind(c(post.sx[draw,]), c(post.sy[draw,])),
                        traploc = traploc, sex = NA, Msexsigma = Msexsigma,
                        logsigmam = post.logsigma[draw], logsigmaf = NA)
      }
      
      loglik = logLfn.piscr.i.waic(logitphi = post.logitphi[draw], logpimat = logpimat, 
                                   z = c(post.z[draw,]), yi00 = yi00,
                                   n0mat = n0mat, n0vec = n0vec, J = J)
      return(loglik)
      
    }}
  
  if(model == 2 | model == 4){
    
    simi.waic = function(draw){
      library(matrixStats)
      right2d.star = right2d[order(post.L[draw,]),]
      if(model == 2){
        
        logpmat = logpmatfn(logitp0 = post.logitp0[draw], s.left = cbind(c(post.sx[draw,]), c(post.sy[draw,])),
                       traploc = traploc, sex = post.sex[draw,], Msexsigma = Msexsigma,
                       logsigmam = post.logsigmam[draw], logsigmaf = post.logsigmaf[draw])
      }
      if(model == 4){
        
        logpmat = logpmatfn(logitp0 = post.logitp0[draw], s.left = cbind(c(post.sx[draw,]), c(post.sy[draw,])),
                       traploc = traploc, sex = NA, Msexsigma = Msexsigma,
                       logsigmam = post.logsigma[draw], logsigmaf = NA)
      }
      
      loglik = logLfn.rpiscr.i.waic(logpmat = logpmat, left2d = left2d,
                                    right2d = right2d.star,
                                    z = c(post.z[draw, ]), J = J)
      return(loglik)
      
    }}
  
#===============================================================================
# Function to compute log-likelihood at each MCMC iteration of every data point
#===============================================================================
 
cps=detectCores(); sfInit(parallel=TRUE, cpus=cps)
sfExportAll(); sfClusterSetupRNG()

start.time=Sys.time()
val.list.u0zs=sfClusterApplyLB(1:(ndraws - burnin),simi.u0zs)
val.list.waic=sfClusterApplyLB(1:(ndraws - burnin),simi.waic)
end.time=Sys.time()
time.taken = end.time - start.time; print(time.taken)
sfStop()

loglik.chain.u0zs=c(unlist(val.list.u0zs))
write.csv(loglik.chain.u0zs, file = paste0(folderpath, '/markovchain.loglikelihood.u0zs.txt', sep = ''), quote = F,row.names = F)

loglik.chain.waic = t(matrix(c(unlist(val.list.waic)), M , tot.length)) # tot.length x M
write.csv(loglik.chain.waic, file = paste0(folderpath, '/markovchain.loglik.waic.txt', sep = ''), quote = F,row.names = F)

} # end of for(model in 1:4)
