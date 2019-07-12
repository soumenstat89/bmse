RPISCRfn = function(piscrobj, 
                    ndraws = 20, burnin = 10, Msexsigma = 1, thin = 1,
                    M = 400, bignum = 10^200, folderpath = NA,
                    nloopL = 50, n.update = 20, dd = 10,
                    batchsize = 1000, mindelta = 0.01, 
                    sigmaOfProposal.logitp0 = 0.08,
                    sigmaOfProposal.logsigmam = 0.02, 
                    sigmaOfProposal.logsigmaf = 0.02,
                    sigmaOfProposal.L = 4, 
                    sigmaOfProposal.s = rep(3, 400),
                    SimstudyIndicator = 1)
{    
  if(Msexsigma == 1){cat('Model = 2', '\n', sep = '')}
  if(Msexsigma == 0){cat('Model = 4', '\n', sep = '')}
  
    xlim = c(unlist(piscrobj$xlim))
    ylim = c(unlist(piscrobj$ylim))
    buffer = piscrobj$buffer
   
    tdf = piscrobj$tdf 
    K = dim(tdf)[1]
    J = dim(tdf)[2] - 3
    dimnames(tdf)[[2]] = c('TrapID', 'Easting', 'Northing', 1:J)
    traploc = as.matrix(tdf[,c('Easting', 'Northing')] )
    mask = tdf[, 4:(J+3)]
    
    #------------------------------------------------------  
    #--------------------------
    # Capture-recapture data
    #--------------------------
    left_edf = piscrobj$edf1;  right_edf = piscrobj$edf2
    if(is.null(left_edf)){ warning("Zero number of captured individuals by detector 1")}
    if(is.null(right_edf)){ warning("Zero number of captured individuals by detector 2")}
    stopifnot(!is.null(left_edf), !is.null(right_edf))
      
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
    numl = numl.obs = dim(left)[1]
    numr = numr.obs = dim(right)[1]
    #------------------------------------
    
    tt = sum(left_edf[,5] == 12)
    if(tt > 0) {IDfixed = left_edf[tt,2]}
    if(tt == 0) {IDfixed = 0}
    
    if (IDfixed == 0) known = 'none'
    if (IDfixed > 0){
      known = 'some'
      if(IDfixed == numl) known = 'ALL'
    }
    
   
    if(Msexsigma == 1)
    {
      sex.data = piscrobj$sex.data
      sexl.data = na.omit(sex.data[,'sex1']) # numl x 1
      sexr.data = na.omit(sex.data[,'sex2']) # numr x 1
      sexl = rep(1, numl)
      sexr = rep(1, numr)
      if(sum(sexl.data == 'Female') > 0){ sexl[sexl.data == 'Female'] = 0} # numl x 1
      if(sum(sexr.data == 'Female') > 0){ sexr[sexr.data == 'Female'] = 0} # numr x 1
   } 
    
    #------------------------------------------------------  
    
    # Also re-label data sets if needed so that the left side always has more individuals
    if (known!="ALL" & numl < numr)
    { 
      a = left; left = right; right = a
      b = numl; numl = numr; numr = b
      if(Msexsigma == 1){ bb = sexl; sexl = sexr; sexr = bb } 
    }
    
    # Using the mask to make some false-positives zero
    for(t in 1:J){ for(k in 1:K){ 
      if(mask[k, t] == 0) {left[1:numl, k, t] = 0; right[1:numr, k, t] = 0}
    }}
    
    left = abind(left, array(0, dim = c( M - numl, K, J)), along = 1)
    right = abind(right, array(0, dim = c( M - numr, K, J)), along = 1)
    
    left2d = apply(left, c(1,2), sum) # no. of times individual i got captured at trap k on LEFT side
    right2d = apply(right, c(1,2), sum) # no. of times individual i got captured at trap k on RIGHT side
    
    #----------------------------------------------
    
    # Initial values
    logitp0 = logit(0.002)
    psi = 0.5
    
    if(Msexsigma == 1)
    {
      logsigmam = log(0.6) # Initialising logsigmaf
      logsigmaf = log(0.5) # Initialising logsigmaf
      theta = 0.6  # Initialising theta

      numlm = sum(sexl == 1, na.rm = T); numlf = sum(sexl == 0, na.rm = T)
      numrm = sum(sexr == 1, na.rm = T); numrf = sum(sexr == 0, na.rm = T)
      
      sexl1 = c(); sexr1 = c()
      if(numlm < numrm) sexl1 = rep(1, numrm - numlm) # sexl = c(sexl, rep(1, numrm - numlm))
      if(numlm > numrm) sexr1 = rep(1, numlm - numrm) # sexr = c(sexr, rep(1, numlm - numrm))
      if(numlf < numrf) sexl1 = c(sexl1, rep(0, numrf - numlf)) # (nnn - numl) x 1 # sexl = c(sexl, rep(0, numrf - numlf))
      if(numlf > numrf) sexr1 = c(sexr1, rep(0, numlf - numrf)) # (nnn - numr) x 1 # sexr = c(sexr, rep(0, numlf - numrf))
      
      sexl = c(sexl, rep(NA, M - numl)) # M x 1
      sexr = c(sexr, rep(NA, M - numr)) # M x 1
      # The guys whose sex are not known or missing
      missing.sex.guys = is.na(sexl) 
      
      nnn = max(numlm, numrm) + max(numlf, numrf) 
      if(numl < nnn){ sexl[(numl+1) : nnn] = sample(sexl1, nnn - numl, replace = F)} # rearranging the new entries in sexl
      if(numr < nnn){ sexr[(numr+1) : nnn] = sample(sexr1, nnn - numr, replace = F)} # rearranging the new entries in sexl
      
      sex.allocation = rbinom(M - nnn, 1, theta)
      sexl[(nnn+1) : M] = sexr[(nnn+1) : M] = sex.allocation # length M and have equal number of males and females
      
      # Latent variable vector for gender
      sex = sexl # M x 1
   } # end of if(Msexsigma == 1)
    
    if(Msexsigma == 0)
    {
      logsigmaf = NA; 
      logsigmam = log(0.6)
      sex = sexl = sexr = NULL
    } 
    #######################################################
    
    
    
    #######################################################
    # 'leftrightmatch' gives initial values for L vector of length M
     if (known != 'ALL')
    {
      L = Lvec.fn(left, right, numl, numr, sexl, sexr, traploc, IDfixed, nloopL, Msexsigma)
    }
    yli00 = rowSums(left) # M x 1 # giving total no. of captures for each LSI
    yri00 = rowSums(right) # M x 1 # giving total no. of captures for each RSI
    
    if (known != 'ALL') {yri00.star = yri00[order(L)];  right.star = right[order(L),,]; right2d.star = right2d[order(L),]}
    if (known == 'ALL') {yri00.star = yri00;            right.star = right; right2d.star = right2d}
    n0mat = apply(left + right.star, c(1,2), function(a){sum(a > 0)}) # trapwise count of captures (in at least one side) for each of the M individual
    n0vec = rowSums(n0mat) 
    y2d = left2d + right2d.star # M x 1 # giving total no. of captures for each captured individual
    yi00 = rowSums(left2d + right2d.star)
    zero.guys = (yi00 == 0) 
    numcap = sum(yi00 > 0) 
    
    # Initializing z
    z = c(rep(1, numcap), rbinom(M - numcap, 1, psi)) 
    
    #==============================
    #Initialising activity centers
    
    s.left = cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
    
    for (i in c(1:M)[yi00>0])
    {
      cap.loc = matrix(c(traploc[y2d[i,] > 0,]), ncol = 2, byrow = F) # get the traps where ith LSI got captured
      if(dim(cap.loc)[1] == 1) s.left[i,] = cap.loc
      if(dim(cap.loc)[1] > 1) s.left[i,] = apply(cap.loc, 2, weighted.mean, w = y2d[i,][y2d[i,] > 0]) 
    }
   
    #==============================
    
    logpmat = logpmatfn(logitp0, s.left, traploc, sex, logsigmam, logsigmaf, Msexsigma)
    llik = logLfn.rpiscr(logpmat, left2d, right2d.star, z, J)
    
    ts = format(Sys.time(), '%d%m%y_%H%M%S')
    
    if(!is.na(folderpath)) # if folderpath is given 
    {
      if(Msexsigma == 1){ folderName = paste(folderpath, '/Bilateral_MCMC_M2_ndraws', ndraws, '_', ts, sep = '')}
      if(Msexsigma == 0){ folderName = paste(folderpath, '/Bilateral_MCMC_M4_ndraws', ndraws, '_', ts, sep = '')}
      dir.create(path = folderName)
    }
    if(is.na(folderpath)) # if folderpath is not given 
    {
      if(Msexsigma == 1){ folderName = paste('Bilateral_MCMC_M2_ndraws', ndraws, '_', ts, sep = '')}
      if(Msexsigma == 0){ folderName = paste('Bilateral_MCMC_M4_ndraws', ndraws, '_', ts, sep = '')}
      dir.create(path = folderName)
    }
    
    continueMCMC = TRUE
    batch = 0
    delta = mindelta
    
    # Setting counters
    naccept.logitp0 = naccept.logsigmam = naccept.logsigmaf = naccept.L = 0
    naccept.s = rep(0, M)
    
    draw = drawp0 = drawsigmam = drawsigmaf = drawL = drawsex = drawact = 0
    
    cat('Begin MCMC sampling:', '\n', '\n')
    
    fname1 = paste(folderName, '/sigma_proposals.txt', sep = '')
    cat(c(paste('sigmaOfProposal.logitp0 = ', sigmaOfProposal.logitp0, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    if(Msexsigma == 1)
    {
      cat(c(paste('sigmaOfProposal.logsigmam = ', sigmaOfProposal.logsigmam, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
      cat(c(paste('sigmaOfProposal.logsigmaf = ', sigmaOfProposal.logsigmaf, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    }
    if(Msexsigma == 0)
    {
      cat(c(paste('sigmaOfProposal.logsigma = ', sigmaOfProposal.logsigmam, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    }
    cat(c(paste('sigmaOfProposal.L = ', sigmaOfProposal.L, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    cat(c(paste('sigmaOfProposal.s = rep(', sigmaOfProposal.s[1], ',M)', sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    
    #--------------------------------------------------------
    cat(c('loglikelihood'), sep = ',', file = paste(folderName, '/markovchain.loglikelihood.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.loglikelihood.txt', sep = ''), append = TRUE)

    if(Msexsigma == 1)
    {
      cat(c('N', 'N.Male', 'p0', 'sigmam', 'sigmaf', 'psi', 'theta', 'numcap'), sep = ',', file = paste(folderName, '/markovchain.txt', sep = ''), append = TRUE)
      cat('\n', file = paste(folderName, '/markovchain.txt', sep = ''), append = TRUE)
      
      cat(c(paste('sex', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.sex.txt', sep = ''), append = TRUE)
      cat('\n', file = paste(folderName, '/markovchain.sex.txt', sep = ''), append = TRUE)
    }
    
    if(Msexsigma == 0)
    {
      cat(c('N', 'p0', 'sigma', 'psi', 'numcap'), sep = ',', file = paste(folderName, '/markovchain.txt', sep = ''), append = TRUE)
      cat('\n', file = paste(folderName, '/markovchain.txt', sep = ''), append = TRUE)
    }
    
    cat(c(paste('z', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.z.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.z.txt', sep = ''), append = TRUE)
    
    cat(c(paste('sx', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.sx.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.sx.txt', sep = ''), append = TRUE)
    
    cat(c(paste('sy', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.sy.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.sy.txt', sep = ''), append = TRUE)
    
    if (known != 'ALL') 
    {
      cat(c(paste('L', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.L.txt', sep = ''), append = TRUE)
      cat('\n', file = paste(folderName, '/markovchain.L.txt', sep = ''), append = TRUE)
    }
    
    start.time = Sys.time()
    while (continueMCMC) {   draw = draw + 1;
    if (draw%%(batchsize/5) == 0)  cat('..... drawing sample #', draw, '\n')
    
    # update the increment/decrement for adaptive Metropolis - Hastings samplers
    if (draw%%batchsize == 0) 
    {
      if (1/sqrt(batch) < mindelta)  delta = 1/sqrt(batch)
      batch = batch + 1
    }
    
    #-----------------
    # update p0 (Uniform prior over (0,1) for p0)
    #-----------------
    
    drawp0 = drawp0 +1
    
    logitp0.cand = rnorm(1, logitp0, sigmaOfProposal.logitp0) 
    logpmat.cand = logpmatfn(logitp0.cand, s.left, traploc, sex, logsigmam, logsigmaf, Msexsigma)
    loglik.cand.p0 = logLfn.rpiscr(logpmat.cand, left2d, right2d.star, z, J)
    loglik.curr.p0 = llik
    
    lognum = loglik.cand.p0   + log(exp(logitp0.cand)/((1+exp(logitp0.cand))^2)) # + dnorm(logitp0, logitp0.cand, sigmaOfProposal.logitp0, log = T))  ## log-prior + log-proposal
    logden = loglik.curr.p0   + log(exp(logitp0)/((1+exp(logitp0))^2)) # + dnorm(logitp0.cand, logitp0, sigmaOfProposal.logitp0, log = T))  ## log-prior + log-proposal
    
    if (logden == -Inf)
    {
      logitp0 = logitp0.cand
      logpmat = logpmat.cand
      loglik.curr.p0 = loglik.cand.p0
      naccept.logitp0 = naccept.logitp0 + 1
    }
    if (logden != -Inf)
    {
      if (runif(1,0,1) <= exp(lognum - logden)) 
      {
        logitp0 = logitp0.cand
        logpmat = logpmat.cand
        loglik.curr.p0 = loglik.cand.p0
        naccept.logitp0 = naccept.logitp0 + 1
      }
    }
    
    if (draw%%batchsize == 0)
    {
      SigmaDiff = ifelse(naccept.logitp0 > 0.44*batchsize, exp(2*delta), exp(-2*delta))
      if(draw <= burnin){ sigmaOfProposal.logitp0 = sigmaOfProposal.logitp0 * SigmaDiff}
      cat("proposal sd of logitp0 = ", sigmaOfProposal.logitp0, ' ')
      cat("naccept.logitp0 = ", naccept.logitp0, '\n')
      naccept.logitp0 = 0   # reset counter for next batch
    }
    
    #-----------------
    # update sigmam ( Uniform prior over (0,dd) for sigmam, dd being maximum stretch of the state space )
    #-----------------
    
    drawsigmam = drawsigmam +1
    
    logsigmam.cand = rnorm(1, logsigmam, sigmaOfProposal.logsigmam)
    logpmat.cand = logpmatfn(logitp0, s.left, traploc, sex, logsigmam.cand, logsigmaf, Msexsigma)
    loglik.cand.sigmam = logLfn.rpiscr(logpmat.cand, left2d, right2d.star, z, J)
    loglik.curr.sigmam = loglik.curr.p0
    
    lognum = loglik.cand.sigmam + log(exp(logsigmam.cand)/dd) 
    logden = loglik.curr.sigmam + log(exp(logsigmam)/dd) 
    
    if (logden == -Inf)
    {
      logsigmam = logsigmam.cand
      logpmat = logpmat.cand
      loglik.curr.sigmam = loglik.cand.sigmam
      naccept.logsigmam = naccept.logsigmam + 1
    }
    if (logden != -Inf)
    {
      if (runif(1,0,1) <= exp(lognum - logden)) 
      {
        logsigmam = logsigmam.cand
        logpmat = logpmat.cand
        loglik.curr.sigmam = loglik.cand.sigmam
        naccept.logsigmam = naccept.logsigmam + 1
      }
    }
    if (draw%%batchsize == 0)
    {
      SigmaDiff = ifelse(naccept.logsigmam > 0.44*batchsize, exp(2*delta), exp(-2*delta))
      if(draw <= burnin){ sigmaOfProposal.logsigmam = sigmaOfProposal.logsigmam * SigmaDiff}
      if(Msexsigma == 1){ 
        cat("proposal sd of logsigmam = ", sigmaOfProposal.logsigmam, ' ')
        cat("naccept.logsigmam = ", naccept.logsigmam, '\n')}
      if(Msexsigma == 0){ 
        cat("proposal sd of logsigma = ", sigmaOfProposal.logsigmam, ' ')
        cat("naccept.logsigma = ", naccept.logsigmam, '\n')}
      naccept.logsigmam = 0   # reset counter for next batch
    }
    
    
    if(Msexsigma == 1){ 
      
      #-----------------
      # update sigmaf ( Uniform prior over (0,dd) for sigmaf, dd being maximum stretch of the state space )
      #-----------------
      
      drawsigmaf = drawsigmaf +1
      
      logsigmaf.cand = rnorm(1, logsigmaf, sigmaOfProposal.logsigmaf)
      logpmat.cand = logpmatfn(logitp0, s.left, traploc, sex, logsigmam, logsigmaf.cand, Msexsigma)
      loglik.cand.sigmaf = logLfn.rpiscr(logpmat.cand, left2d, right2d.star, z, J)
      loglik.curr.sigmaf = loglik.curr.sigmam
      
      lognum = loglik.cand.sigmaf  + log(exp(logsigmaf.cand)/dd) 
      logden = loglik.curr.sigmaf  + log(exp(logsigmaf)/dd) 
      
      if (logden == -Inf) 
      {
        logsigmaf = logsigmaf.cand
        logpmat = logpmat.cand
        loglik.curr.sigmaf = loglik.cand.sigmaf
        naccept.logsigmaf = naccept.logsigmaf + 1
      }
      if (logden != -Inf)
      {
        if (runif(1,0,1) <= exp(lognum - logden)) 
        {
          logsigmaf = logsigmaf.cand
          logpmat = logpmat.cand
          loglik.curr.sigmaf = loglik.cand.sigmaf
          naccept.logsigmaf = naccept.logsigmaf + 1
        }
      } 
      if (draw%%batchsize == 0)
      {
        SigmaDiff = ifelse(naccept.logsigmaf > 0.44*batchsize, exp(2*delta), exp(-2*delta))
        if(draw <= burnin){ sigmaOfProposal.logsigmaf = sigmaOfProposal.logsigmaf * SigmaDiff}
        cat("proposal sd of logsigmaf = ", sigmaOfProposal.logsigmaf, ' ')
        cat("naccept.logsigmaf = ", naccept.logsigmaf, '\n')
        naccept.logsigmaf = 0   # reset counter for next batch
      }
    } # end of if(Msexsigma == 1)   
    
    #-----------------
    #  update L (Uniform prior over permutation space of 1:M)
    #-----------------
    
    if (known != "ALL")
    {
      loglik.curr.L = ifelse(Msexsigma == 1, loglik.curr.sigmaf, loglik.curr.sigmam)
      drawL = drawL +1
      indx1 = (IDfixed + 1):M
      
      if( IDfixed < numr) rightset1  = order(L)[z == 1 & order(L) >= (IDfixed+1) & order(L) <= numr]
      if( IDfixed >= numr) rightset1  = c()
      
      rightset2 = order(L)[z == 1 & order(L) > numr] # real right individuals wh went uncaptured
      if(length(rightset2) >1) rightset2 = sample(rightset2, min(length(rightset2), n.update), replace = F)
      rightset4 = unique(c(rightset1, rightset2))
      rightset = sort(sample(rightset4, min(length(rightset4), n.update), replace = F)) 
     
      for (r.guy1 in rightset)
      {
        
        l.swap.out = L[r.guy1]
        
        if(Msexsigma == 1){ possible.L = c(1:M)[z == 1 & c(1:M) > IDfixed & sex == sex[l.swap.out]]} # Real LSI's with same gender as l.swap.out among the non-IDfixed
        if(Msexsigma == 0){ possible.L = c(1:M)[z == 1 & c(1:M) > IDfixed]} # Real LSI's among the non-IDfixed
        
        dv =  (s.left[l.swap.out,1] - s.left[,1]) ^ 2 + (s.left[l.swap.out,2] - s.left[,2]) ^ 2 
        wt.possible.L = exp( - dv / sigmaOfProposal.L ^ 2)[possible.L]
        
        if (length(possible.L) > 1)  l.swap.in =  sample( possible.L, 1, replace = F, prob = wt.possible.L)
        if (length(possible.L) == 1) { naccept.L = naccept.L + 1; next}  # It means, there does not exist any LSI in the neighbourhood of l.swap.out except itself.
        if (length(possible.L) == 0) { next } # this case will never happen since l.swap.out is present there at the centre of the circle dv < 5
        if (l.swap.in == l.swap.out) {naccept.L = naccept.L + 1; next} # saves computation time in for loop
      
        ##  Which right encounter history is currently associated with left guy s.swap.in?
        r.guy2 = c(1:M)[L == l.swap.in]
        L.cand = L; L.cand[r.guy1] = l.swap.in; L.cand[r.guy2] = l.swap.out
        dv.back = (s.left[l.swap.in,1] - s.left[,1]) ^ 2 + (s.left[l.swap.in,2] - s.left[,2]) ^ 2  
        wt.possible.back.L = exp( - dv.back / sigmaOfProposal.L ^ 2)[possible.L]
        
        jump.prob.L =  wt.possible.L[which(possible.L == l.swap.in, arr.ind = T)] / sum(wt.possible.L) #  q(state.curr, state.cand)
        jump.prob.back.L =  wt.possible.back.L[which(possible.L == l.swap.out, arr.ind = T)] / sum(wt.possible.back.L) #  q(state.cand, state.curr)
        right.star.cand = right[order(L.cand),,]
        right2d.star.cand = apply(right.star.cand, c(1,2), sum)
        
        loglik.cand.L = logLfn.rpiscr(logpmat, left2d, right2d.star.cand, z, J)
        
        lognum = loglik.cand.L + log(jump.prob.back.L)
        logden = loglik.curr.L + log(jump.prob.L)
        
        if (logden == -Inf) 
        {
          L = L.cand 
          loglik.curr.L = loglik.cand.L
          right.star = right.star.cand
          right2d.star = right2d.star.cand
          naccept.L = naccept.L + 1
        }
        if (logden != -Inf)
        {
          if (runif(1,0,1) <= exp(lognum - logden)) 
          {
            L = L.cand 
            loglik.curr.L = loglik.cand.L
            right.star = right.star.cand
            right2d.star = right2d.star.cand
            naccept.L = naccept.L + 1
          } 
        }
        
      } # end of for (r.guy1 in rightset)
      
      if (draw%%batchsize == 0)
      {
        SigmaDiff = ifelse(naccept.L > 0.44*batchsize, exp(2*delta), exp(-2*delta))
        if(draw <= burnin){ sigmaOfProposal.L= sigmaOfProposal.L * SigmaDiff}
        cat(paste("proposal sd of L = ", sep = ""), sigmaOfProposal.L, ' ')
        cat(paste("naccept.L = ", sep = ""), naccept.L, '\n')
        naccept.L = 0   # reset counter for next batch
      } 
      
    } # end of if (known!='ALL')
    
    #======================================
    yi00 = rowSums(left2d + right2d.star) 
    numcap = sum( yi00 > 0) 
    zero.guys = (yi00 == 0) 
    
    ncprob = (1 - exp(logpmat)) ^ (2*J)  # M x K
    prob0 = exp(rowSums(log(ncprob)))  # M x 1
    
    fc = prob0*psi / (prob0*psi + 1 - psi)
    z[zero.guys] = rbinom(sum(zero.guys), 1, fc[zero.guys])
    z[!zero.guys] = 1
    
    #-----------------
    # update psi (Uniform prior over (0,1))
    #-----------------
    psi = rbeta(1, 1 + sum(z), 1 +  (M - sum(z) ) )  # count IDfixed guys
    #-----------------
    # update the activity centers (Uniform prior over state space S)
    #-----------------
    drawact = drawact + 1
    
    s.left.cand = matrix(c(rnorm(M, s.left[, 1], sigmaOfProposal.s), rnorm(M, s.left[, 2], sigmaOfProposal.s)), nrow = M, ncol = 2 , byrow = F)
    inbox = s.left.cand[,1] < xlim[2] & s.left.cand[,1] > xlim[1] & s.left.cand[,2] < ylim[2] & s.left.cand[,2] > ylim[1]
    logpmat.cand = logpmatfn(logitp0, s.left.cand, traploc, sex, logsigmam, logsigmaf, Msexsigma)
    
    # simi2 = function(i) 
    for(ind in c(1:M)[z == 1 & inbox])
    {
      part1 = dbinom(left2d[ind,], J, exp(logpmat.cand[ind,]), log = T) # K x 1
      part2 = dbinom(right2d.star[ind,], J, exp(logpmat.cand[ind,]), log = T) # K x 1
      part1[part1 == -Inf] = -bignum; part2[part2 == -Inf] = -bignum; 
      loglik.cand.s = z[ind] * sum(part1 + part2)
      
      part1 = dbinom(left2d[ind,], J, exp(logpmat[ind,]), log = T) # K x 1
      part2 = dbinom(right2d.star[ind,], J, exp(logpmat[ind,]), log = T) # K x 1
      part1[part1 == -Inf] = -bignum; part2[part2 == -Inf] = -bignum; 
      loglik.curr.s = z[ind] * sum(part1 + part2)
      
      lognum = loglik.cand.s 
      logden = loglik.curr.s 
      
      if (logden == -Inf)
      {
        s.left[ind, ] = s.left.cand[ind, ]
        naccept.s[ind] = naccept.s[ind] + 1
      }
      if (logden != -Inf)
      {
        if (runif(1,0,1) <= exp(lognum - logden)) 
        {
          s.left[ind, ] = s.left.cand[ind, ]
          naccept.s[ind] = naccept.s[ind] + 1
        }
      }
      if (draw%%batchsize == 0) 
      {
        SigmaDiff = ifelse(naccept.s[ind] > 0.44*batchsize, exp(2*delta), exp(-2*delta))
        if(draw <= burnin){ sigmaOfProposal.s[ind] = sigmaOfProposal.s[ind] * SigmaDiff}
        naccept.s[ind] = 0   # reset counter for next batch
      } 
      
    } # end of   for(i in c(1:M)[z == 1 & inbox])
    
    #-----------------
    #  update sex (Bernoulli(theta) prior)
    #-----------------
    if(Msexsigma == 1) 
    {  
      drawsex = drawsex +1
      
      D1 = e2dist(s.left, traploc)
      logpmat.m = log(expit(logitp0)) - D1 * D1 / (2 * (exp(logsigmam) ^ 2)) # M x K # Division is always column wise
      part1.m = rowSums(dbinom(left2d, J, exp(logpmat.m), log = T)) # M x 1
      part2.m = rowSums(dbinom(right2d.star, J, exp(logpmat.m), log = T)) # M x 1
      part1.m[part1.m == -Inf] = -bignum; part2.m[part2.m == -Inf] = -bignum; # .Machine$double.xmax
      AAm = exp(part1.m + part2.m)
      
      logpmat.f = log(expit(logitp0)) - D1 * D1 / (2 * (exp(logsigmaf) ^ 2)) # M x K # Division is always column wise
      part1.f = rowSums(dbinom(left2d, J, exp(logpmat.f), log = T)) # M x 1
      part2.f = rowSums(dbinom(right2d.star, J, exp(logpmat.f), log = T)) # M x 1
      part1.f[part1.f == -Inf] = -bignum; part2.f[part2.f == -Inf] = -bignum; # .Machine$double.xmax
      AAf = exp(part1.f + part2.f)
      
      AA = theta * AAm + (1-theta) * AAf # M x 1
      theta0 = theta * AAm / AA # M x 1
      mis = (z == 1) & missing.sex.guys
      sex[mis] = rbinom(sum(mis), 1, theta0[mis])
      
      #-----------------
      # update theta  (Uniform prior over state space (0,1))
      #-----------------
      theta = rbeta( 1, 1 + sum(z*sex), 1 + sum(z*(1 - sex)) )
      
    } # end of if(Msexsigma == 1) 
    
    logpmat = logpmatfn(logitp0, s.left, traploc, sex, logsigmam, logsigmaf, Msexsigma)
    llik = logLfn.rpiscr(logpmat, left2d, right2d.star, z, J)
 
    cat(llik, sep = ',', file = paste(folderName, '/markovchain.loglikelihood.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.loglikelihood.txt', sep = ''), append = TRUE)
    
    if(Msexsigma == 1)
    {
      cat(c(sum(z), sum(sex*z), expit(logitp0), 
            exp(logsigmam), #*scale,
            exp(logsigmaf), #*scale
            psi, theta, numcap
      ), sep = ',', file = paste(folderName, '/markovchain.txt', sep = ''), append = TRUE)
      cat('\n', file = paste(folderName, '/markovchain.txt', sep = ''), append = TRUE)
      
      cat(c(sex), sep = ',', file = paste(folderName, '/markovchain.sex.txt', sep = ''), append = TRUE) # note that the first numl elements of sexl are observed hence kept fixed in the markov chain.
      cat('\n', file = paste(folderName, '/markovchain.sex.txt', sep = ''), append = TRUE)
    }  # end of if(Msexsigma == 1)
    
    if(Msexsigma == 0)
    {
      cat(c(sum(z), expit(logitp0), 
            exp(logsigmam), #*scale,
            psi, numcap
      ), sep = ',', file = paste(folderName, '/markovchain.txt', sep = ''), append = TRUE)
      cat('\n', file = paste(folderName, '/markovchain.txt', sep = ''), append = TRUE)
      
    } # end of if(Msexsigma == 0)
    
    cat(c(z), sep = ',', file = paste(folderName, '/markovchain.z.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.z.txt', sep = ''), append = TRUE)
    
    cat(c(s.left[,1]), sep = ',', file = paste(folderName, '/markovchain.sx.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.sx.txt', sep = ''), append = TRUE)
    
    cat(c(s.left[,2]), sep = ',', file = paste(folderName, '/markovchain.sy.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.sy.txt', sep = ''), append = TRUE)
    
    if (known != 'ALL')
    {
      cat(c(L), sep = ',', file = paste(folderName, '/markovchain.L.txt', sep = ''), append = TRUE)
      cat('\n', file = paste(folderName, '/markovchain.L.txt', sep = ''), append = TRUE)
    }
    #-------------------------------------------  
    if (draw == ndraws) 
    {
      cat('Completed ', ndraws, ' draws of MCMC algorithm', '\n')
      numOfDraws = 0      
      if (numOfDraws == 0) continueMCMC = FALSE
      
      cat(c(paste('\n', 'Completed ', ndraws, ' draws of MCMC algorithm', sep = ''), '\n \n'), sep = ',', file = fname1, append = TRUE)
      cat(c(paste('sigmaOfProposal.logitp0 = ', sigmaOfProposal.logitp0, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
      if(Msexsigma == 1)
      {
        cat(c(paste('sigmaOfProposal.logsigmam = ', sigmaOfProposal.logsigmam, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
        cat(c(paste('sigmaOfProposal.logsigmaf = ', sigmaOfProposal.logsigmaf, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
      }
      if(Msexsigma == 0)
      {
        cat(c(paste('sigmaOfProposal.logsigma = ', sigmaOfProposal.logsigmam, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
      }
      cat(c(paste('sigmaOfProposal.L = ', sigmaOfProposal.L, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
      cat(c(paste('sigmaOfProposal.s = rep(', sigmaOfProposal.s[1], ',M)', sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
      
    } # end of if (draw == ndraws) 
    #-------------------------------------------   
    }  # end of loop for MCMC draws
    
    cat('MCMC sampling is completed!', '\n', '\n')
    end.time = Sys.time()
    (time.taken = end.time - start.time); print(time.taken)
    #==================================================================================
    
    # Compute estimates of population abundance for a specific region within S
    
    # .... compute Bayes estimate of posterior mean N and estimate of posterior probabilities
    
    post = read.csv(paste(folderName, '/markovchain.txt', sep = ''), sep = ',', header = T)
    seq0 = seq(burnin+1, ndraws, by = thin)
    post = post[seq0,]
    post = post[, -ncol(post)]
    geweke_diag = geweke.diag(post, frac1 = 0.2, frac2 = 0.4)
    print(geweke_diag)
    suppressWarnings(write.table(round(unlist(geweke_diag), 3), sep = '           ', 
                                 col.names = F, append = F,
                                 file = paste(folderName, '/geweke.diag.txt', sep = '')))
    
    N.chain = post[, 'N']
    
    Nvalues = 0:M
    probN = rep(0, (M + 1))
    for (i in 1:(M + 1)) probN[i] = length(N.chain[N.chain == (i - 1)])/((ndraws - burnin)/thin)
    
    post.mode.N = Nvalues[probN == max(probN)][1] 
    Bayes.Nmean = mean(N.chain, na.rm = T)
    Bayes.Nvar = mean(N.chain ^ 2, na.rm = T) - (mean(N.chain, na.rm = T)) ^ 2
    ind1 = cumsum(probN) >= 0.025 & cumsum(probN) <= 0.975
    Bayes.Nlower = quantile(N.chain, 0.025, na.rm = T, type = 1) 
    Bayes.Nupper = quantile(N.chain, 0.975, na.rm = T, type = 1) 
    # # # .... plot histograms of empirical Bayes and Bayes estimates
    fname5 = paste(folderName, '/mcmc plots of N.jpeg', sep = '')
    ylimits = c(0, max(probN))
    jpeg(fname5, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
    plot(Nvalues[1:M], probN[1:M], type = 'h', ylim = ylimits, xlab = 'N', ylab = 'probability' )
    dev.off()
    
    #=====================================#
    # Print results of Bayesian analysis  #
    #=====================================#
    out = as.matrix(c(Bayes.Nmean, sqrt(Bayes.Nvar), Bayes.Nlower, post.mode.N, Bayes.Nupper))
    dimnames(out) = list(c('Mean.N', 'SE.N', '2.5%', 'post.mode', '97.5%'),c('HB Estimates of N'))
    prob.quantiles = c(0.025, 0.5, 0.975)  # for credible limits
    post.stats = cbind(apply(post,2,mean, na.rm = T), apply(post,2,sd, na.rm = T), 
                       t(apply(post, 2, quantile, probs = prob.quantiles, na.rm = T,type = 1)))
    
    if(SimstudyIndicator == 1){
      
      Msexsigma.given = piscrobj$Msexsigma.given 
      N.given = piscrobj$N.given
      phi.given = piscrobj$phi.given
      omega0.given = piscrobj$omega0.given
      if(Msexsigma.given == 1)
      {
        N.Male.given = piscrobj$N.Male.given
        sigmam.given = piscrobj$sigma.given[1]
        sigmaf.given = piscrobj$sigma.given[2]
      } 
      if(Msexsigma.given == 0)
      {
        N.Male.given = sigmaf.given = NA
        sigmam.given = piscrobj$sigma.given[1]
      } 
      
    if(Msexsigma == 1)
    { 
      par_true = c(N.given, N.Male.given, omega0.given, 
                   sigmam.given, sigmaf.given, N.given / M, N.Male.given / N.given)
      post.stats2 = t(cbind(apply(rbind(par_true, post), 2, sim.stats, hpdprob = 0.95)))
    }
    
    if(Msexsigma == 0)
    {
      par_true = c(N.given, omega0.given, 
                   sigmam.given, N.given / M)
      post.stats2 = t(cbind(apply(rbind(par_true, post), 2, sim.stats, hpdprob = 0.95)))
    }
      dimnames(post.stats2)[2] = list(c('bias', 'relative_bias_percentage', 'RMSE', 'HPD_lower', 'HPD_upper'))
      
    }
   
    prob.names = paste(as.character(100*prob.quantiles), '%', sep = '')
    dimnames(post.stats)[2] = list(c('Mean.Chain', 'SD.Chain', prob.names))
    
    #============================================================================================================
    #=============================#
    #     Monte carlo errors      #
    #=============================#
    mcse.mean.vec = mcse.sd.vec = rep(0, dim(post)[2])
    mcse.lq.vec = mcse.uq.vec = mcse.med.vec = rep(0, dim(post)[2])
    for (i in 1:dim(post)[2])
    {
      postvec = post[,i]
      mcse.mean.vec[i] = mcse(postvec[!is.na(postvec)], mean)$se
      mcse.sd.vec[i] = mcse(postvec[!is.na(postvec)], sd)$se
      mcse.med.vec[i] = mcse(postvec[!is.na(postvec)], medQuantile)$se
      mcse.lq.vec[i] = mcse(postvec[!is.na(postvec)], lowerQuantile)$se
      mcse.uq.vec[i] = mcse(postvec[!is.na(postvec)], upperQuantile)$se
    }
    mcse.mean.vec = unlist(mcse.mean.vec)
    mcse.sd.vec = unlist(mcse.sd.vec)
    mcse.med.vec = unlist(mcse.med.vec)
    mcse.lq.vec = unlist(mcse.lq.vec)
    mcse.uq.vec = unlist(mcse.uq.vec)
    mcse.mat = cbind(mcse.mean.vec, mcse.sd.vec, mcse.lq.vec, mcse.med.vec, mcse.uq.vec)
    
    HB_estimates = sim_estimates = MCSE_estimates = c('', '', '', '', '')
    dim.names1 = c('Mean.Chain', 'SD.Chain', prob.names)
    dim.names2 = c('bias', 'relative_bias_percentage', 'RMSE', 'HPD_lower', 'HPD_upper')
    dimnames(mcse.mat) = dimnames(post.stats)
    
    if(SimstudyIndicator == 1){
 
    information = as.data.frame(c(M, N.given, N.Male.given, IDfixed, numl.obs, numr.obs, known, phi.given, omega0.given, 
                                  sigmam.given, sigmaf.given, K, J, ndraws, burnin, thin, buffer,
                                xlim, ylim, Msexsigma, Msexsigma.given))
    dimnames(information) = list(c('M', 'N.given', 'N.Male.given', 'IDfixed', 'num1', 'num2', 'known', 'phi.given', 'omega0.given', 
                                 'sigmam.given', 'sigmaf.given', 'ntrap', 'nocc', 'ndraws', 'burnin', 'thin', 'buffer',
                                 'xlim.1', 'xlim.2', 'ylim.1', 'ylim.2', 'Msexsigma', 'Msexsigma.given'),
                                 c('info'))
    write.csv(t(information), file = paste(folderName, '/info1.csv', sep = ''), quote = F,row.names = T)
    
    out.final = rbind(t(out), HB_estimates, dim.names1, post.stats, sim_estimates, dim.names2, post.stats2, MCSE_estimates, dim.names1, mcse.mat)
    cat('Parameter estimates', '\n', sep = '')
    print(post.stats)
    print(post.stats2)
    }
    
    if(SimstudyIndicator == 0){
      information = as.data.frame(c(M, IDfixed, numl.obs, numr.obs, known, K, J, ndraws, burnin, thin, buffer,
                                    xlim, ylim, Msexsigma))
      dimnames(information) = list(c('M', 'IDfixed', 'num1', 'num2', 'known', 'ntrap', 'nocc', 'ndraws', 'burnin', 'thin', 'buffer',
                                     'xlim.1', 'xlim.2', 'ylim.1', 'ylim.2', 'Msexsigma'),
                                   c('info'))
      write.csv(t(information), file = paste(folderName, '/info1.csv', sep = ''), quote = F,row.names = T)
      
      out.final = rbind(t(out), HB_estimates, dim.names1, post.stats, sim_estimates, MCSE_estimates, dim.names1, mcse.mat)
      cat('Parameter estimates', '\n', sep = '')
      print(post.stats)
    }
    cat('\n', 'Monte Carlo Error', '\n', sep = '')
    print(mcse.mat)
    
    write.csv(rbind(cbind(out.final, matrix('',nrow(out.final), nrow(information) - ncol(out.final))), 
                    dimnames(information)[[1]], t(information)), 
              file = paste(folderName, '/EstimatesOfDerivedParam.csv', sep = ''), quote = F,row.names = T)
    
    
    #============================================================================================================
    #=====================#
    #     Traceplots      #
    #=====================#
    post = read.csv(paste(folderName, '/markovchain.txt', sep = ''), sep = ',', header = T)
    post.mcmc = as.mcmc(post)
    
    fname6 = paste(folderName, '/traceplots_N.jpeg', sep = '')
    jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
    traceplot(post.mcmc[, 'N'], xlab = 'Iterations', ylab = 'N', main = 'Traceplot of N',
              cex.main = 2, cex.lab = 2, cex.axis = 2)
    dev.off()
    
    fname6 = paste(folderName, '/traceplots_psi.jpeg', sep = '')
    jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
    traceplot(post.mcmc[, 'psi'], xlab = 'Iterations', ylab = 'psi', main = 'Traceplot of psi',
              cex.main = 2, cex.lab = 2, cex.axis = 2)
    dev.off()
    
    fname6 = paste(folderName, '/traceplots_p0.jpeg', sep = '')
    jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
    traceplot(post.mcmc[, 'p0'], xlab = 'Iterations', ylab = 'p0', main = 'Traceplot of p0',
              cex.main = 2, cex.lab = 2, cex.axis = 2)
    dev.off()
    
    if(Msexsigma == 0)
    {  
      fname6 = paste(folderName, '/traceplots_sigma.jpeg', sep = '')
      jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
      traceplot(post.mcmc[, 'sigma'], xlab = 'Iterations', ylab = 'sigma.male', main = 'Traceplot of sigma.male',
                cex.main = 2, cex.lab = 2, cex.axis = 2)
      dev.off()
    }
    
    if(Msexsigma == 1)
    {
      fname6 = paste(folderName, '/traceplots_N.Male.jpeg', sep = '')
      jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
      traceplot(post.mcmc[, 'N.Male'], xlab = 'Iterations', ylab = 'N.Males', main = 'Traceplot of N.Males',
                cex.main = 2, cex.lab = 2, cex.axis = 2)
      dev.off()
      
      fname6 = paste(folderName, '/traceplots_theta.jpeg', sep = '')
      jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
      traceplot(post.mcmc[, 'theta'], xlab = 'Iterations', ylab = 'theta', main = 'Traceplot of theta',
                cex.main = 2, cex.lab = 2, cex.axis = 2)
      dev.off()
      
      fname6 = paste(folderName, '/traceplots_sigmam.jpeg', sep = '')
      jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
      traceplot(post.mcmc[, 'sigmam'], xlab = 'Iterations', ylab = 'sigma.male', main = 'Traceplot of sigma.male',
                cex.main = 2, cex.lab = 2, cex.axis = 2)
      dev.off()
      
      fname6 = paste(folderName, '/traceplots_sigmaf.jpeg', sep = '')
      jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
      traceplot(post.mcmc[, 'sigmaf'], xlab = 'Iterations', ylab = 'sigma.female', main = 'Traceplot of sigma.female',
                cex.main = 2, cex.lab = 2, cex.axis = 2)
      dev.off()
    }
    #============================================================================================================
    #=====================#
    #    Scatter plots    #
    #=====================#
    post.mcmc = as.mcmc(post)
    
    post = read.csv(paste(folderName, '/markovchain.txt', sep = ''), sep = ',', header = T)
    seq0 = seq(burnin+1, ndraws, by = thin)
    post = post[seq0,]
    seq1 = seq(1, dim(post)[1], by = 10) 
    N1 = post[seq1,'N']; N1 = N1[!is.na(N1)]
    psi1 = post[seq1,'psi']; psi1 = psi1[!is.na(psi1)]
    p01 = post[seq1,'p0']; p01 = p01[!is.na(p01)]
    
    if(Msexsigma == 0)
    {
      sigma1 = post[seq1,'sigma']; sigma1 = sigma1[!is.na(sigma1)]
      df = data.frame(N1, psi1, p01, sigma1) 
      
      fname6 = paste(folderName, '/Scatter Plots_', IDfixed,known, '.jpeg', sep = '')
      jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
      pairs( ~ N1 + psi1 + p01 + sigma1, data = df, main='')
      dev.off()
    }
    
    if(Msexsigma == 1)
    {
      N.Male1 = post[seq1, 'N.Male']; N.Male1 = N.Male1[!is.na(N.Male1)]
      theta1 = post[seq1, 'theta']; theta1 = theta1[!is.na(theta1)]
      sigmam1 = post[seq1,'sigmam']; sigmam1 = sigmam1[!is.na(sigmam1)]
      sigmaf1 = post[seq1, 'sigmaf']; sigmaf1 = sigmaf1[!is.na(sigmaf1)]
      df = data.frame(N1, psi1, N.Male1, theta1, p01, sigmam1, sigmaf1) 
      
      fname6 = paste(folderName, '/Scatter Plots_', IDfixed,known, '.jpeg', sep = '')
      jpeg(fname6, width = 1000, height = 1000, units = 'px', pointsize = 12, quality = 100)
      pairs( ~ N1 + psi1 + N.Male1 + theta1 + p01 + sigmam1 + sigmaf1, data = df, main='')
      dev.off()
    }
    
    #====================================================
write.csv(piscrobj$edf1, file = paste0(folderName, '/encounter_data_file1.csv', sep = ''), quote = F,row.names = F)
write.csv(piscrobj$edf2, file = paste0(folderName, '/encounter_data_file2.csv', sep = ''), quote = F,row.names = F)
write.csv(piscrobj$tdf, file = paste0(folderName, '/trap_deployment_file.csv', sep = ''), quote = F,row.names = F)
if(Msexsigma.given == 1){write.csv(piscrobj$sex.data, file = paste0(folderName, '/sex_captured.csv', sep = ''), quote = F,row.names = F)}
   #====================================================
    save.image(paste(folderName, '/savingRimage.RData', sep = ''))

} # end of RPISCRfn
