
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
   
      ##======================
      ## Log-likelihood chain
      ##======================
 
      loglik.chain.waic = read.csv(paste0(folderpath, '/markovchain.loglik.waic.txt', sep = ''), sep = ',', header = T) # tot.length x M

      if(model == 1 | model  == 2)
      {
        logfactor.sex = log(post.theta^(post.z*post.sex)) + log((1 - post.theta)^(1 - post.z*post.sex)) # tot.length x M
        loglik.chain.waic = loglik.chain.waic + logfactor.sex # tot.length x M
      }

      loglik.chain.dic = rowSums(loglik.chain.waic) # tot.length x 1
     
      indxvec = unlist(read.csv(paste(folderpath, '/iter_indxinit_indxvec.csv', sep = ''), sep = ',', header = T, row.names = 1))
      
      iter = c(indxvec[1])
      initindx = c(indxvec[2])
      if(length(indxvec) == 2){indxvec = c()} 
      if(length(indxvec) > 2){indxvec = c(indxvec[3:length(indxvec)])} # length(indxvec) = iter
      
      #============================================================
      
      if(iter%%2 == 1){
        if(iter == 1){ indxx1 = indxx2 = initindx}
        if(iter > 1){ indxx1 = indxvec[iter-2]; indxx2 = indxvec[iter-1]}} # length(indxvec) = iter
      
      if(iter%%2 == 0){
        if(iter == 2){ indxx1 = indxvec[iter-1]; indxx2 = initindx}
        if(iter > 2){ indxx1 = indxvec[iter-1]; indxx2 = indxvec[iter-2]}}
      
      psi.star = post.psi[indxx1]; 
      if(model == 1 | model  == 3) { phi.star = post.phi[indxx1]; omega0.star = post.omega0[indxx1]; }
      if(model == 2 | model  == 4) { p0.star = post.p0[indxx1]; }
      if(model == 1 | model  == 2) { theta.star = post.theta[indxx1]; sigmam.star = post.sigmam[indxx1]; sigmaf.star = post.sigmaf[indxx1]}
      if(model == 3 | model  == 4) { sigma.star = post.sigma[indxx1]}
      if(model == 1 | model  == 2) { sex.star = post.sex[indxx2,]}
      z.star = post.z[indxx2,]; sx.star = post.sx[indxx2,]; sy.star = post.sy[indxx2,]; L.star = post.L[indxx2,]
      if(model == 1 | model == 3){
       
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
            loglik.pt = loglik.pt0 + logfactor.sex 
          }
          
          if(model == 3)
          { # Dey et al. (2019) Partial ID WITHOUT GENDER
            logpimat = logpimatfn(logitbasepar = logit(omega0.star), s.left = cbind(sx.star, sy.star), traploc = traploc, sex = NULL,
                            logsigmam = log(sigma.star), logsigmaf = NA, Msexsigma = Msexsigma)
            loglik.pt0 = logLfn.piscr(logitphi = logit(phi.star), logpimat = logpimat, z = z.star,
                                      yi00 = yi00, n0mat = n0mat, n0vec = n0vec, J = J)
            loglik.pt = loglik.pt0 
          }
         
        } 
      
      if(model == 2 | model == 4){
        
          right2d.star = right2d[order(L.star),]
          if(model == 2)
          { # ROYLE (2015) Partial ID WITH GENDER
            logpmat = logpmatfn(logitp0 = logit(p0.star), s.left = cbind(sx.star, sy.star), traploc = traploc, sex = sex.star,
                           logsigmam = log(sigmam.star), logsigmaf = log(sigmaf.star), Msexsigma = Msexsigma)
            loglik.pt0 = logLfn.rpiscr.dic(logpmat = logpmat, left2d = left2d, right2d = right2d.star, z = z.star, J = J)
            logfactor.sex = sum(log(theta.star^(z.star*sex.star))) + sum(log((1 - theta.star)^(1 - z.star*sex.star)))
            loglik.pt = loglik.pt0 + logfactor.sex 
          }
          
          if(model == 4)
          { # ROYLE (2015) Partial ID WITHOUT GENDER
            logpmat = logpmatfn(logitp0 = logit(p0.star), s.left = cbind(sx.star, sy.star), traploc = traploc, sex = NULL,
                           logsigmam = log(sigma.star), logsigmaf = NA, Msexsigma = Msexsigma)
            loglik.pt0 = logLfn.rpiscr.dic(logpmat = logpmat, left2d = left2d, right2d = right2d.star, z = z.star, J = J)
            loglik.pt = loglik.pt0 
          }
         
        }
   
      #=======================================
      # DIC
      #=======================================
      
      val.argmax = dicfn(loglik.chain.dic, loglik.pt)
      dic.val.argmax = val.argmax[1:2]
      p_dic = val.argmax[3:4]
      
      out = cbind(c(dic.val.argmax), c(p_dic))
      
      dimnames(out) = list(c('DIC1', 'DIC2'), c('COL1', 'COL2'))
     
      write.csv(out, file = paste(folderpath, '/dic.csv', sep = ''), quote = F, row.names = T)
      
    } # end of for(model in 1:4)

#====================================================
save.image(paste0(folderpath, '/savingRimage_DIC.RData', sep = ''))
time.taken = Sys.time() - start.time; print(time.taken)
#=================================================
