
rm(list = ls())

codefolder = getwd()
ts = format(Sys.time(), "%d%m%y_%H%M%S")
folderName = paste(codefolder, '/BayesianModelSelection_', ts, sep = '')
dir.create(path = folderName)
setwd(folderName)

# Loading the required functions
source(paste(codefolder, '/BMSE.utility.functions.R', sep = ''))
source(paste(codefolder, '/BMSE.PISCRfn.R', sep = ''))
source(paste(codefolder, '/BMSE.Royle.PISCRfn.R', sep = ''))

# Simulating data set
piscrobj = sim.partial.data(N = 100, N.Male = 40, phi = 0.5, omega0 = 0.005,
sigma = c(0.3, 0.15), Msexsigma = 1, nrow_trap = 16, ncol_trap = 10,
nocc = 50, xlim = c(0,5), ylim = c(0,7), buffer = 1)

# Running MCMC
for(i in c(1,0)){

# MCMC as per Dey et al. model
PISCRfn(piscrobj, Msexsigma = i, ndraws = 20, burnin = 10, thin  = 1,  M = 400,
        bignum = 10^200, folderpath = NA, nloopL = 50, n.update = 20, dd = 10,
        batchsize = 1000, mindelta = 0.01, sigmaOfProposal.logitphi = 0.05,
        sigmaOfProposal.logitomega0 = 0.08, sigmaOfProposal.logsigmam = 0.02,
        sigmaOfProposal.logsigmaf = 0.02, sigmaOfProposal.L = 4,
        sigmaOfProposal.s = rep(3, 400), SimstudyIndicator = 1)

# MCMC as per Royle (2015) model
RPISCRfn(piscrobj, Msexsigma = i, ndraws = 20, burnin = 10, thin  = 1, M = 400,  
         bignum = 10^200, folderpath = NA, nloopL = 50, n.update = 20, dd = 10, 
         batchsize = 1000, mindelta = 0.01, sigmaOfProposal.logitp0 = 0.08,
         sigmaOfProposal.logsigmam = 0.02, sigmaOfProposal.logsigmaf = 0.02,
         sigmaOfProposal.L = 4, sigmaOfProposal.s = rep(3, 400), SimstudyIndicator = 1)
}

# Computing integrated likelihood and log-pointwise likelihood
source(paste(codefolder, '/BMSE.intlik.u0zs.waic.R', sep = ''))

# Estimating log-marginal likelihood with integrate likelihood approach 
source(paste(codefolder, '/GDBF.IL.R', sep = ''))

# Estimating log-marginal likelihood with MAP approach 
source(paste(codefolder, '/GDBF.MAP.R', sep = ''))

# Computing DIC
source(paste(codefolder, '/BMSE.DIC.R', sep = ''))

# Computing WAIC
source(paste(codefolder, '/BMSE.WAIC.R', sep = ''))

# Computing posterior predictive loss
source(paste(codefolder, '/BMSE.PPL.R', sep = ''))


# Obtaining table of values for the obtained model selection criteria
folderNames1 = c()
for(model in 1:4) { folderNames1 = c(folderNames1, paste(getwd(), '/', fff0[grep(paste("_M", model, "_", sep = ""), fff0)], '', sep = ''))}

for(model in 1:4){
 
  cat('Model = ', model, '\n', sep = '')
  folderpath = folderNames1[model]
 
  margmap = read.csv(paste(folderpath, '/marginal.MAP.csv', sep = ''), sep = ',', header = T, row.names = 1)
  margil = read.csv(paste(folderpath, '/marginal.IL.csv', sep = ''), sep = ',', header = T, row.names = 1)
  dic = read.csv(paste(folderpath, '/dic.csv', sep = ''), sep = ',', header = T, row.names = 1)
  waic = read.csv(paste(folderpath, '/waic.csv', sep = ''), sep = ',', header = T, row.names = 1)
  ppl = read.csv(paste(folderpath, '/ppl.csv', sep = ''), sep = ',', header = T, row.names = 1)
  msc = rbind(margmap, margil, dic, waic, ppl)
  tuning_den_names.MAP = c('MAP.GD.Normal',  paste('MAP.GD.t', c(10, 100, 500, 1000, 10000), sep = ''),
                           paste('MAP.GD.Trunc_normal_conf', c(0.9,0.95, 0.99), sep = ''))
  tuning_den_names.IL = c('IL.GD.Normal',  paste('IL.GD.t', c(10, 100, 500, 1000, 10000), sep = ''),
                          paste('IL.GD.Trunc_normal_conf', c(0.9,0.95, 0.99), sep = ''))
  
  dimnames(msc) = list(c(tuning_den_names.MAP,'HM',tuning_den_names.IL, 'DIC1', 'DIC2', 'WAIC1', 'WAIC2', 'WAIC3', 'PPL'),
                       c('COL1', 'COL2'))
  write.csv(msc, file = paste(folderpath, '/all_model_selection_criteria.csv', sep = ''), quote = F, row.names = T)
  print(msc)
} 
