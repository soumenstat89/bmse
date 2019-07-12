#=================================================================
# Function to compute MAP (maximum a posteriori)
#=================================================================
MAPfn = function(data, loglikfn, logpriorfn = NULL, par1.chain, par2.chain, 
                 loglik.chain, logprior.chain = NULL, burnin = NA){
  # Please ensure par1.chain and par2.chain are either a vector or a matrix or a data frame
  # loglikfn = It is a function to compute loglikelihood. 
  #            The arguments of loglikfn are data, par1, par2 (in that order), 
  #            e.g., loglikfn(data, par1, par2)
  # logpriorfn = It is a function to compute prior density.
  #              The arguments of logpriorfn are par1, par2 (in that order),
  #              e.g., logpriorfn(par1, par2)
  # loglik.chain = it is a vector of loglikelihod values for each MCMC iterations
  # logprior.chain = it is vector of logprior density values for each MCMC iterations
  # par1.chain = a matrix of MCMC chain for par1 with chain of 
  #             different parameters in different columns
  # par2.chain = a matrix of MCMC chain for par2 with chain of 
  #             different parameters in different columns
  # data = It should be a single R object (e.g., list, matrix,
  #        data.frame etc.) which takes all the non-paramter elements
  #        (i.e., capture history, mask data, trap deployment data).
  #        and to be supplied as argument of the likelihood function
  # burnin = burnin value of the MCMC chain. 
  #          In case the chains are already truncated, please provide
  #          the value of `burnin' as 0 (zero).
  library(doMC)
  ndraws = length(loglik.chain)
  if(is.na(burnin)){burnin = ndraws/2}
  tot.length = ndraws - burnin
  par1.chain = as.matrix(par1.chain); par1.chain = par1.chain[(burnin + 1):ndraws,]
  par2.chain = as.matrix(par2.chain); par2.chain = par2.chain[(burnin + 1):ndraws,]
  loglik.chain = loglik.chain[(burnin+1):ndraws]
  if(!is.null(logprior.chain)){
    logprior.chain = logprior.chain[((burnin+1):ndraws)]
    logLP.chain = logLP.chain.init = loglik.chain + logprior.chain
  }
  if(is.null(logprior.chain)){
    logLP.chain = logLP.chain.init = loglik.chain
  }
  newmax.logLP.chain = initmax.logLP = max(logLP.chain)
  indxnewmax = initindx = c(1:tot.length)[logLP.chain == max(logLP.chain)][1] # which(loglik.chain.dic==max(loglik.chain.dic), arr.ind = T)
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
       cat(' Iteration (MAP) = ', iter, '\n', sep = '')
       cps=detectCores() 
       registerDoMC(cps - 1)
       val.list <- foreach(draw=1:tot.length) %dopar% {                       
          if(iter%%2 == 1){
            par1.star = par1.chain[draw,]
            par2.star = par2.chain[indxx,]
          }
          if(iter%%2 == 0){ 
            par1.star = par1.chain[indxx,]
            par2.star = par2.chain[draw,]
          }
          if(!is.null(logpriorfn)){
            loglikprior = loglikfn(data, par1.star, par2.star) +
                         logpriorfn(par1.star, par2.star)
          }
          if(is.null(logpriorfn)){
            loglikprior = loglikfn(data, par1.star, par2.star) 
          }
          loglikprior
        } 
        newlogLP.chain=c(unlist(val.list))
        newmax.logLP.chain = max(newlogLP.chain)
        indxnewmax = c(1:tot.length)[newlogLP.chain == newmax.logLP.chain][1] 
        indxvec = c(indxvec, indxnewmax)
        maxLPvec = c(maxLPvec, newmax.logLP.chain)
        condition = newmax.logLP.chain > max.logLP.chain
  } # end of while(condition & iter < tot.length) loop
  cat('Iterations have finished. \n', sep = '')
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
  #============================================================
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
  #============================================================
  out = list(par1.map = par1.chain[indxx1, ], par2.map = par2.chain[indxx2, ],
             iterdata = t(iterdata), loglikprior = c(initmax.logLP, maxLPvec))
  return(out)
} # end of function MAPfn
      
      
    
