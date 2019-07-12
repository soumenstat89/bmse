
# rm(list = ls())
options(digits = 8)

start.time = Sys.time()
ts = format(Sys.time(), "%d%m%y_%H%M%S")

iter.star = rep(NA, 4)
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
      state.space.area = (xlim[2] - xlim[1])*(ylim[2] - ylim[1]) # = 5 * 7 = 35 = data$state.space.area
      R = sqrt((xlim[2]-xlim[1])^2 + (ylim[2]-ylim[1])^2) # R = Maximum stretch of the state space
      
      #======================
      # Log-prior chain 
      #======================
      # Prior (Transformed) 
      # <[Priors are transformed because tuning density uses transformed parameters]>
      logprior.z = rowSums(log(post.psi^post.z) + log((1 - post.psi)^(1 - post.z))) # tot.length x 1
      logprior.s = log(1/state.space.area) # scalar
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
        logprior.logsigmam + logprior.logsigmaf + logprior.s + logprior.L  }
      if(model == 2){ logprior.chain = logprior.logitpsi + logprior.logittheta + logprior.logitp0 + 
        logprior.logsigmam + logprior.logsigmaf + logprior.s + logprior.L }
      if(model == 3){ logprior.chain = logprior.logitpsi + logprior.logitphi + logprior.logitomega0 + 
        logprior.logsigma + logprior.s + logprior.L }
      if(model == 4){ logprior.chain = logprior.logitpsi + logprior.logitp0 + 
        logprior.logsigma + logprior.s + logprior.L }
      
      if(model == 1){ post.chains = cbind(post.logitpsi, post.logittheta, post.logitphi, post.logitomega0, post.logsigmam, post.logsigmaf)}
      if(model == 2){ post.chains = cbind(post.logitpsi, post.logittheta, post.logitp0, post.logsigmam, post.logsigmaf) }
      if(model == 3){ post.chains = cbind(post.logitpsi, post.logitphi, post.logitomega0, post.logsigma) }
      if(model == 4){ post.chains = cbind(post.logitpsi, post.logitp0, post.logsigma) }
      ###############################################################
      
      loglik.chain =  unlist(read.csv(paste0(folderpath, '/markovchain.loglikelihood.txt', sep = ''), sep = ',', header = T))[(burnin + 1):ndraws]
      if(model == 1 | model  == 2)
      {
        logfactor.sex = log(post.theta^(post.z*post.sex)) + log((1 - post.theta)^(1 - post.z*post.sex)) # tot.length x M
        loglik.chain = loglik.chain + rowSums(logfactor.sex)  # tot.length x 1
      }
      
      logLP.chain = logLP.chain.init = loglik.chain + logprior.chain + logprior.z 
      newmax.logLP.chain = initmax.logLP = max(logLP.chain)
      indxnewmax = initindx = c(1:tot.length)[logLP.chain == max(logLP.chain)][1] # which(loglik.chain.dic==max(loglik.chain.dic), arr.ind = T)
    
      #============================
      # To compute MAP estimate
      #============================
    
      indxvec = c()      
      maxLPvec = c()
      iter = 0
      condition = TRUE  
      cat('MAP iteration begins. \n', sep = '')
      while(condition & iter < (tot.length^2 - 1))
      {
       
       iter = iter + 1
       indxx = indxnewmax
       
       max.logLP.chain = newmax.logLP.chain
       cat(' Model = ', model, ' Iteration = ', iter, '\n', sep = '')

      if(model == 1 | model == 3){
        
        simi.map = function(draw){
         
          if(iter%%2 == 1){
            psi.star = post.psi[draw]; phi.star = post.phi[draw]
            if(model == 1 | model  == 3) { omega0.star = post.omega0[draw]}
            if(model == 2 | model  == 4) { p0.star = post.p0[draw]}
            if(model == 1 | model  == 2) { theta.star = post.theta[draw]; sigmam.star = post.sigmam[draw]; sigmaf.star = post.sigmaf[draw]}
            if(model == 3 | model  == 4) { sigma.star = post.sigma[draw]}
            if(model == 1 | model  == 2) { sex.star = post.sex[indxx,]}
            z.star = post.z[indxx,]; sx.star = post.sx[indxx,]; sy.star = post.sy[indxx,]; L.star = post.L[indxx,]}
          
          if(iter%%2 == 0){
            psi.star = post.psi[indxx]; phi.star = post.phi[indxx]
            if(model == 1 | model  == 3) { omega0.star = post.omega0[indxx]}
            if(model == 2 | model  == 4) { p0.star = post.p0[indxx]}
            if(model == 1 | model  == 2) { theta.star = post.theta[indxx]; sigmam.star = post.sigmam[indxx]; sigmaf.star = post.sigmaf[indxx]}
            if(model == 3 | model  == 4) { sigma.star = post.sigma[indxx]}
            if(model == 1 | model  == 2) { sex.star = post.sex[draw,]}
            z.star = post.z[draw,]; sx.star = post.sx[draw,]; sy.star = post.sy[draw,]; L.star = post.L[draw,]}
          
          right.star = right[order(L.star),,]
          right2d.star = right2d[order(L.star),]
          n0mat = apply(left + right.star, c(1,2), function(a){sum(a > 0)}) # trapwise count of recorded entries for each of the M individual
          n0vec = rowSums(n0mat) # giving total count of recorded entries for each of the M individual
          yli00 = rowSums(left2d); yri00.star = rowSums(right2d.star); yi00 = yli00 + yri00.star 
          if(model == 1)
          { # Dey et al. (2019) Partial ID WITH GENDER
            logpimat = logpimatfn(logitbasepar = logit(omega0.star), s.left =  cbind(sx.star, sy.star), traploc = traploc, sex = sex.star,
                            logsigmam = log(sigmam.star), logsigmaf = log(sigmaf.star), Msexsigma = Msexsigma)
            loglik.pt0 = logLfn.piscr(logitphi = logit(phi.star), logpimat = logpimat, z = z.star,
                                      yi00 = yi00, n0mat = n0mat, n0vec = n0vec, J = J)
            logfactor.sex = sum(log(theta.star^(z.star*sex.star))) + sum(log((1 - theta.star)^(1 - z.star*sex.star)))
            logfactor.z = sum(log(psi.star^z.star) + log((1 - psi.star)^(1 - z.star))) 
            loglik.pt = loglik.pt0 + logfactor.sex + logfactor.z
          }
          
          if(model == 3)
          { # Dey et al. (2019) Partial ID WITHOUT GENDER
            logpimat = logpimatfn(logitbasepar = logit(omega0.star), s.left = cbind(sx.star, sy.star), traploc = traploc, sex = NULL,
                            logsigmam = log(sigma.star), logsigmaf = NA, Msexsigma = Msexsigma)
            loglik.pt0 = logLfn.piscr(logitphi = logit(phi.star), logpimat = logpimat, z = z.star,
                                     yi00 = yi00, n0mat = n0mat, n0vec = n0vec, J = J)
            logfactor.z = sum(log(psi.star^z.star) + log((1 - psi.star)^(1 - z.star))) 
            loglik.pt = loglik.pt0 + logfactor.z
          }
          return(loglik.pt)
          
        }}
      
      if(model == 2 | model == 4){
        simi.map = function(draw){
          if(iter%%2 == 1){
            psi.star = post.psi[draw]; 
            if(model == 1 | model  == 3) { omega0.star = post.omega0[draw]}
            if(model == 2 | model  == 4) { p0.star = post.p0[draw]}
            if(model == 1 | model  == 2) { theta.star = post.theta[draw]; sigmam.star = post.sigmam[draw]; sigmaf.star = post.sigmaf[draw]}
            if(model == 3 | model  == 4) { sigma.star = post.sigma[draw]}
            if(model == 1 | model  == 2) { sex.star = post.sex[indxx,]}
            z.star = post.z[indxx,]; sx.star = post.sx[indxx,]; sy.star = post.sy[indxx,]; L.star = post.L[indxx,]}
          
          if(iter%%2 == 0){
            psi.star = post.psi[indxx]; 
            if(model == 1 | model  == 3) { omega0.star = post.omega0[indxx]}
            if(model == 2 | model  == 4) { p0.star = post.p0[indxx]}
            if(model == 1 | model  == 2) { theta.star = post.theta[indxx]; sigmam.star = post.sigmam[indxx]; sigmaf.star = post.sigmaf[indxx]}
            if(model == 3 | model  == 4) { sigma.star = post.sigma[indxx]}
            if(model == 1 | model  == 2) { theta.star = post.sex[draw,]}
            z.star = post.z[draw,]; sx.star = post.sx[draw,]; sy.star = post.sy[draw,]; L.star = post.L[draw,]}
          right2d.star = right2d[order(L.star),]
          if(model == 2)
          { # ROYLE (2015) Partial ID WITH GENDER
            logpmat = logpmatfn(logitp0 = logit(p0.star), s.left = cbind(sx.star, sy.star), traploc = traploc, sex = sex.star,
                           logsigmam = log(sigmam.star), logsigmaf = log(sigmaf.star), Msexsigma = Msexsigma)
            loglik.pt0 = logLfn.rpiscr.dic(logpmat = logpmat, left2d = left2d, right2d = right2d.star, z = z.star, J = J)
            logfactor.sex = sum(log(theta.star^(z.star*sex.star))) + sum(log((1 - theta.star)^(1 - z.star*sex.star)))
            logfactor.z = sum(log(psi.star^z.star) + log((1 - psi.star)^(1 - z.star))) 
            loglik.pt = loglik.pt0 + logfactor.sex + logfactor.z
          }
          
          if(model == 4)
          { # ROYLE (2015) Partial ID WITHOUT GENDER
            logpmat = logpmatfn(logitp0 = logit(p0.star), s.left = cbind(sx.star, sy.star), traploc = traploc, sex = NULL,
                           logsigmam = log(sigma.star), logsigmaf = NA, Msexsigma = Msexsigma)
            loglik.pt0 = logLfn.rpiscr.dic(logpmat =logpmat, left2d = left2d, right2d = right2d.star, z = z.star, J = J)
            logfactor.z = sum(log(psi.star^z.star) + log((1 - psi.star)^(1 - z.star))) 
            loglik.pt = loglik.pt0 + logfactor.z
          }
          return(loglik.pt)
        }}  
        cps=detectCores(); sfInit(parallel=TRUE, cpus=  cps) 
        sfExportAll(); sfClusterSetupRNG(); restoreResults <- TRUE
        start.time=Sys.time()
        val.list=sfClusterApplySR(1:tot.length,simi.map,perUpdate=10,restore=restoreResults)
        end.time=Sys.time()
        time.taken = end.time - start.time; print(time.taken)
        sfStop()
        
        loglik.chain.bf=c(unlist(val.list))
        
          newlogLP.chain = loglik.chain.bf + logprior.chain
          newmax.logLP.chain = max(newlogLP.chain)
          indxnewmax = c(1:tot.length)[newlogLP.chain == newmax.logLP.chain][1] 
          if(iter%%2 == 1){outlogLP.chain = newlogLP.chain}
          
        indxvec = c(indxvec, indxnewmax)
        maxLPvec = c(maxLPvec, newmax.logLP.chain)
        condition = newmax.logLP.chain > max.logLP.chain
        
     } # end of # while(condition)
     cat('MAP iteration finished. \n', sep = '')
     
     if(length(indxvec) > 0){
     iterdata = as.data.frame(c(iter, initindx, indxvec))
     dimnames(iterdata) = list(c('niter.map', 'mcmc.init.indx', 
                                 paste('mcmcindx', 1:length(indxvec), sep = '')),
                               c('iterdata'))
     }
     
     if(length(indxvec) == 0){
       iterdata = as.data.frame(c(iter, initindx))
       dimnames(iterdata) = list(c('niter.map', 'mcmc.init.indx'),
                                 c('iterdata'))
     }
     
     write.csv(t(iterdata), file = paste(folderpath, '/iter_indxinit_indxvec.csv', sep = ''), quote = F,row.names = T)
     cat(c(initmax.logLP, maxLPvec), file = paste0(folderpath, '/initmax.logLP_maxLPvec.txt', sep = ''), sep = ',')
     #============================================================
     if(iter < (tot.length^2 - 1)){
       if(iter%%2 == 1){
         if(iter == 1){ indxx1 = indxx2 = initindx}
         if(iter > 1){  indxx1 = indxvec[iter-2]; indxx2 = indxvec[iter-1]}
       } 
       
       if(iter%%2 == 0){
         if(iter == 2){ indxx1 = indxvec[iter-1]; indxx2 =  initindx}
         if(iter > 2){  indxx1 = indxvec[iter-1]; indxx2 = indxvec[iter-2]}
       }
     }
     
     if(iter == (tot.length^2 - 1)){
       if(iter%%2 == 1){
         if(condition == F){ indxx1 = indxvec[iter-2]; indxx2 = indxvec[iter-1]}
         if(condition == T){ indxx1 = indxvec[iter]; indxx2 = indxvec[iter-1]}
       } 
       
       if(iter%%2 == 0){
         if(condition == F){ indxx1 = indxvec[iter-1]; indxx2 = indxvec[iter-2]}
         if(condition == T){ indxx1 = indxvec[iter-1]; indxx2 = indxvec[iter]}
       }
     }
     
     if(model == 1 | model == 3){
       simi.map = function(draw){
        
           psi.star = post.psi[draw]; phi.star = post.phi[draw]
           if(model == 1 | model  == 3) { omega0.star = post.omega0[draw]}
           if(model == 2 | model  == 4) { p0.star = post.p0[draw]}
           if(model == 1 | model  == 2) { theta.star = post.theta[draw]; sigmam.star = post.sigmam[draw]; sigmaf.star = post.sigmaf[draw]}
           if(model == 3 | model  == 4) { sigma.star = post.sigma[draw]}
           if(model == 1 | model  == 2) { sex.star = post.sex[indxx2,]}
           z.star = post.z[indxx2,]; sx.star = post.sx[indxx2,]; sy.star = post.sy[indxx2,]; L.star = post.L[indxx2,]
           
         right.star = right[order(L.star),,]
         right2d.star = right2d[order(L.star),]
         n0mat = apply(left + right.star, c(1,2), function(a){sum(a > 0)}) # trapwise count of recorded entries for each of the M individual
         n0vec = rowSums(n0mat) # apply(n0mat, 1, sum) # giving total count of recorded entries for each of the M individual
         yli00 = rowSums(left2d); yri00.star = rowSums(right2d.star); yi00 = yli00 + yri00.star
         if(model == 1)
         { # Dey et al. (2019) Partial ID WITH GENDER
           logpimat = logpimatfn(logitbasepar = logit(omega0.star), s.left =  cbind(sx.star, sy.star), traploc = traploc, sex = sex.star,
                           logsigmam = log(sigmam.star), logsigmaf = log(sigmaf.star), Msexsigma = Msexsigma)
           loglik.pt0 = logLfn.piscr(logitphi = logit(phi.star), logpimat = logpimat, z = z.star,
                                     yi00 = yi00, n0mat = n0mat, n0vec = n0vec, J = J)
           logfactor.sex = sum(log(theta.star^(z.star*sex.star))) + sum(log((1 - theta.star)^(1 - z.star*sex.star)))
            loglik.pt = loglik.pt0 + logfactor.sex 
         }

         if(model == 3)
         { # Dey et al. (2019) Partial ID WITHOUT GENDER
           logpimat = logpimatfn(logitbasepar= logit(omega0.star), s.left = cbind(sx.star, sy.star), traploc = traploc, sex = NULL,
                           logsigmam = log(sigma.star), logsigmaf = NA, Msexsigma = Msexsigma)
           loglik.pt0 = logLfn.piscr(logitphi = logit(phi.star), logpimat = logpimat, z = z.star,
                                     yi00 = yi00, n0mat = n0mat, n0vec = n0vec, J = J)
           loglik.pt = loglik.pt0 
         }
         return(loglik.pt)
         
       }}

     if(model == 2 | model == 4){
       simi.map = function(draw){
         
           psi.star = post.psi[draw]; 
           if(model == 1 | model  == 3) { omega0.star = post.omega0[draw]}
           if(model == 2 | model  == 4) { p0.star = post.p0[draw]}
           if(model == 1 | model  == 2) { theta.star = post.theta[draw]; sigmam.star = post.sigmam[draw]; sigmaf.star = post.sigmaf[draw]}
           if(model == 3 | model  == 4) { sigma.star = post.sigma[draw]}
           if(model == 1 | model  == 2) { sex.star = post.sex[indxx2,]}
           z.star = post.z[indxx2,]; sx.star = post.sx[indxx2,]; sy.star = post.sy[indxx2,]; L.star = post.L[indxx2,]

         right2d.star = right2d[order(L.star),]
         if(model == 2)
         { # ROYLE (2015) Partial ID WITH GENDER
           logpmat = logpmatfn(logitp0 = logit(p0.star), s.left = cbind(sx.star, sy.star), traploc = traploc, sex = sex.star,
                          logsigmam = log(sigmam.star), logsigmaf = log(sigmaf.star), Msexsigma = Msexsigma)
           loglik.pt0 = logLfn.rpiscr.dic(logpmat = logpmat, left2d = left2d, right2d = right2d.star, z = z.star, J = J)
           logfactor.sex = sum(log(theta.star^(z.star*sex.star))) + sum(log((1 - theta.star)^(1 - z.star*sex.star)))
           loglik.pt = loglik.pt0 + logfactor.sex # + logfactor.z
         }

         if(model == 4)
         { # ROYLE (2015) Partial ID WITHOUT GENDER
           logpmat = logpmatfn(logitp0 = logit(p0.star), s.left = cbind(sx.star, sy.star), traploc = traploc, sex = NULL,
                          logsigmam = log(sigma.star), logsigmaf = NA, Msexsigma = Msexsigma)
           loglik.pt0 = logLfn.rpiscr.dic(logpmat = logpmat, left2d = left2d, right2d = right2d.star, z = z.star, J = J)
           loglik.pt = loglik.pt0 
         }
         return(loglik.pt)
         
       }}
     cps=detectCores(); sfInit(parallel=TRUE, cpus= cps)
     sfExportAll(); sfClusterSetupRNG(); restoreResults <- TRUE
     start.time=Sys.time()
     val.list.bf=sfClusterApplySR(1:tot.length,simi.map,perUpdate=10,restore=restoreResults)
     end.time=Sys.time()
     time.taken = end.time - start.time; print(time.taken)
     sfStop()

     loglik.chain.bf=c(unlist(val.list.bf))
     
     if(model == 1){ logprior.chain.bf = logprior.logitpsi + logprior.logittheta + logprior.logitphi + logprior.logitomega0 + 
       logprior.logsigmam + logprior.logsigmaf } 
     if(model == 2){ logprior.chain.bf = logprior.logitpsi + logprior.logittheta + logprior.logitp0 + 
       logprior.logsigmam + logprior.logsigmaf } 
     if(model == 3){ logprior.chain.bf = logprior.logitpsi + logprior.logitphi + logprior.logitomega0 + 
       logprior.logsigma } 
     if(model == 4){ logprior.chain.bf = logprior.logitpsi + logprior.logitp0 + 
       logprior.logsigma } 
     
     logLP.chain = loglik.chain.bf + logprior.chain.bf
     #==================================
     # Parameters for tuning function
     #==================================
      
      mu.est = apply(post.chains, 2, mean, na.rm = T)
      Sigma.est.choice2 = ((dim(post.chains)[1]-1)/dim(post.chains)[1])*cov(post.chains)
      diag(Sigma.est.choice2)[diag(Sigma.est.choice2)< 10^(-5)] = 10^(-5)
      out2 = log.g.chain2 = NULL
      
      #==============================================
      # Gelfand-Dey marginal likelihood estimator
      #==============================================
      log.g.norm.chain = dmvnorm(x = post.chains, mean = mu.est,
                                 sigma = Sigma.est.choice2, log = T) 
      log.g.chain2 = cbind(log.g.chain2, log.g.norm.chain)
      logh.normalg.chain = log.g.norm.chain-logLP.chain 
      logmarglik.normalg = gdmean(logh.normalg.chain)
      out2 = rbind(out2, c(logmarglik.normalg,NA))
      
      degf = c(10, 100, 500, 1000, 10000)
     
      for(i in 1:length(degf)){
        log.g.tg.chain = dmvt(x = post.chains, delta = mu.est, sigma = Sigma.est.choice2,
                              df =degf[i], log = T, type = "shifted")   
        log.g.chain2 = cbind(log.g.chain2, log.g.tg.chain)
        logh.tg.chain = log.g.tg.chain- logLP.chain 
        logmarglik.tg = gdmean(logh.tg.chain)
        out2 = rbind(out2, c(logmarglik.tg,NA))
      } # end of for loop
     
      conflev = c(0.9, 0.95, 0.99)
      q.chisq = qchisq(conflev, df = dim(post.chains)[2], ncp = 0, lower.tail = T)
      scoreval = apply(post.chains, 1, scorefn, mu = mu.est, Sigma = Sigma.est.choice2)
      
      for(lev in 1:length(conflev)){
        indices = scoreval < q.chisq[lev]
        log.g.gewekeg.chain = dmvtruncnorm(x = post.chains[indices,], mu = mu.est,
                                              Sigma = Sigma.est.choice2,
                                              conf.lev = conflev[lev], logg = T) # + logprior.L
        if(length(log.g.gewekeg.chain) == tot.length) {log.g.chain2 = cbind(log.g.chain2, log.g.gewekeg.chain)}
        if(length(log.g.gewekeg.chain) < tot.length) {log.g.chain2 = cbind(log.g.chain2, c(log.g.gewekeg.chain, 
                                                                                           rep(NA, tot.length - length(log.g.gewekeg.chain))))}
        logh.gewekeg.chain = log.g.gewekeg.chain- logLP.chain[indices] #   (loglik.zx0s.chain[indices]+logprior.chain[indices])
        logmarglik.gewekeg = gdmean(logh.gewekeg.chain)
        out2 = rbind(out2, c(logmarglik.gewekeg, NA))
      } # end of for loop
     
      #=================================================
      # Harmonic mean estimator of marginal likelihood
      #=================================================
      # g(mu, L) = pi(mu, L) i.e Takng prior of (mu, L) as the tuning density of (mu, L) 
      loglik.chain =  unlist(read.csv(paste0(folderpath, '/markovchain.loglikelihood.txt', sep = ''), sep = ',', header = T))[(burnin + 1):ndraws]
      if(model == 1 | model  == 2){
        logfactor.sex = log(post.theta^(post.z[,1:numl]*post.sex[,1:numl]))  + log((1 - post.theta)^(1 - post.z[,1:numl]*post.sex[,1:numl])) # tot.length x numl
        loglik.chain = loglik.chain + rowSums(logfactor.sex)  # tot.length x 1
      }
      logh.chain = - loglik.chain 
      logmarglik.hm = gdmean(logh.chain)
      out2 = rbind(out2, c(logmarglik.hm,NA))
      
      tuning_den_names = c('GD.Normal',  paste('GD.t', c(10, 100, 500, 1000, 10000), sep = ''),
                           paste('GD.Trunc_normal_conf', c(0.9,0.95, 0.99), sep = ''))
      
      dimnames(out2)= list(c(tuning_den_names, 'HM'), c('COL1', 'COL2'))
      dimnames(log.g.chain2)[[2]] = tuning_den_names
      write.csv(out2, file = paste(folderpath, '/marginal.MAP.csv', sep = ''), quote = F, row.names = T)
      write.csv(log.g.chain2, file = paste(folderpath, '/chain.log.g.MAP.txt', sep = ''), quote = F,row.names = F)
      write.csv(logLP.chain, file = paste(folderpath, '/logLik_logPrior.chain.MAP.txt', sep = ''), quote = F, row.names = F)
  
    } # end of for(model in 1:4)
 
#=================================================
save.image(paste0(folderpath, '/savingRimage_GDBF_MAP.RData', sep = ''))
time.taken = Sys.time() - start.time; print(time.taken)
#=================================================

