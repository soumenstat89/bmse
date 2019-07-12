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
  

  #===============================================================
  
simi = function(draw)
{
  if(model == 1){ phi1 = post.phi[draw]; basepar1 = post.omega0[draw]; sigma1 = c(post.sigmam[draw], post.sigmaf[draw]); sex1 = post.sex[draw,]}
  if(model == 2){ phi1 = NULL; basepar1 = post.p0[draw]; sigma1 = c(post.sigmam[draw], post.sigmaf[draw]); sex1 = post.sex[draw,]}
  if(model == 3){ phi1 = post.phi[draw]; basepar1 = post.omega0[draw]; sigma1 = post.sigma[draw]; sex1 = NULL}
  if(model == 4){ phi1 = NULL; basepar1 = post.p0[draw]; sigma1 = post.sigma[draw]; sex1 = NULL}
  
  data = sim.partial.data.ppl(phi = phi1, basepar = basepar1,
             sigma = sigma1, sex = sex1, 
             z = post.z[draw,], L = post.L[draw,], 
             S = cbind(post.sx[draw,], post.sy[draw,]),
             trap.locations = traploc, nocc = J, model = model) 
  return(data)  
}

    obsd.data = c(c(left), c(right)) # 2MKJ x 1 
    data.sim = simi(1) 
    sim.sum.left = data.sim$crdata1; sim.sumsq.left = data.sim$crdata1^2; # M x K x J
    sim.sum.right = data.sim$crdata2; sim.sumsq.right = data.sim$crdata2^2; # M x K x J
    for( draw in 2:tot.length)
    { 
      data.sim = simi(draw)
      sim.sum.left = sim.sum.left + data.sim$crdata1;  sim.sumsq.left = sim.sumsq.left + data.sim$crdata1^2 # M x K x J
      sim.sum.right = sim.sum.right + data.sim$crdata2; sim.sumsq.right = sim.sumsq.right + data.sim$crdata2^2 # M x K x J
    }
    mean.pred.vec = c(c(sim.sum.left), c(sim.sum.right)) / tot.length # 2MKJ x 1
    var.pred.vec = c(c(sim.sumsq.left), c(sim.sumsq.right)) / tot.length - mean.pred.vec^2 # 2MKJ x 1
    ppl.val = sum(c(obsd.data - mean.pred.vec)^2) + sum(var.pred.vec)
    
    out = cbind(c(ppl.val), NA)
    
    dimnames(out) = list(c('PPL'), c('COL1', 'COL2'))
    
    write.csv(out, file = paste(folderpath, '/ppl.csv', sep = ''), quote = F, row.names = T)
   
} # end of for(model in 1:nmodel)
#====================================================
save.image(paste0(folderpath, '/savingRimage_PPL.RData', sep = ''))
time.taken = Sys.time() - start.time; print(time.taken)
#=================================================
