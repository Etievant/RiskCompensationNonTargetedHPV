## -----------------------------------------------------------------------------
## Description: This script illustrates estimation of the NDE when using the 
##              methods proposed in Etievant et al. (Biometrics, 2023) in an 
##              observational setting with risk compensation
##
##              The simulation is displayed in Section 4 in the Main Document of 
##              Etievant (2025)
##
##              Comments:
##              Treatment (i.e., vaccine) assignment has not been randomized
##              Some confounders have not been observed
##              The vaccine has a direct and an indirect effect on the targeted 
##              outcome because of risk compensation. The mediator is unmeasured
## -----------------------------------------------------------------------------

### load packages and files ----------------------------------------------------
source("EstimationFunctions.R")   # functions to be used for the estimation
load("Parameters.RData")          # file with some of the fixed parameters
library(parallel)
library(ggplot2)
library(xtable)
library(openxlsx)

### Notations ------------------------------------------------------------------
#   Y1      infection with the HPV strain targeted by the vaccine (e.g., HPV 16). 
#           We have Y1i = 1 if the ith person is infected with this type, and 0 
#           otherwise. This is the argeted (i.e., primary) outcome we consider
#
#   Y2      contains the Y2^(j), i.e., infection with HPV type j not targeted by 
#           the vaccine. We have Y2i^(j) = 1 if the ith person is infected with 
#           this type, and 0 otherwise
#
#   Y2s     categorical variable indicating the total number of infections by 
#           HPV types not targeted by the vaccine, i.e., Y2s = sum_j Y2^(j). 
#           This is the nonargeted (i.e., secondary) outcome we consider
#
#   T       treatment, that is HPV vaccination. We have Ti = 1 if the ith person
#           has been vaccinated, and 0 otherwise
#
#   A       unobserved confounder e.g. sexual activity pre-vaccination
#
#   W       observed confounder e.g. age and geographical region
#
#   Atilda  unobserved mediator e.g. sexual activity after vaccination

### Varying parameters ---------------------------------------------------------
P_Y1  = rbind(0.14, 0.05, 0.025) # incidence of vaccine targeted HPV type

Aa    = rbind(c(0, 1, 2.5), c(0, 1, 2), c(0, 0.75, 1.5)) # sets of 3 categories 
# for A, that will induce different levels of correlation between Y1 and Y2. We
# Someone who is not sexually active cannot be infected by any of the HPV types 
# (first category is zero)

N     = c(5000, 10000) # observational study sizes

combinations    = expand.grid(n = N, Aa_row = 1:nrow(Aa), P_Y1 = P_Y1)
PARAM           = cbind(combinations$P_Y1, Aa[combinations$Aa_row, ], combinations$n) # parameter values for the different scenarios
colnames(PARAM) = c("P_Y1", "alow", "amedium", "ahigh", "n")

### Set fixed parameters -------------------------------------------------------
wRegion         = c(0,1,2) # an individual can come from 1 of 3 regions 
nRegion         = length(wRegion)
PwRegion        = rep(1/nRegion, nRegion) # P(Wregion = wregion) = 1/3 for each region

nAge            = 13
wAge            = seq(15, 21, length.out = nAge) # 13 age groups
PwAge           = param$PwAge
PwAge           = PwAge / rowSums(PwAge) # P(Wage = wage) for each wage in wAge
# we use slightly different probabilities for each of the k strata region

Pa        = array(NA, c(nRegion, nAge, 3)) # A can take 3 values 
# (alow, amedium and ahigh), with probability depending on the value of W
Pa[,,1]   = param$Palow    # P(A = alow | W)
Pa[,,2]   = param$Pamedium  # P(A = amedium | W)
Pa[,,3]   = param$Pahigh    # P(A = ahigh | W)

# parameters used for the simulation of T given A and W
delta_T   = - 1 / 18
gamma_T   = 1.5
zeta_T    = 1
alpha_T   = -0.91

a                   = Aa[1, ] # Atilda can also take values alow, amedium and 
# ahigh, with probability depending on the values of W, A and T
Patildalow          = array(0, c(nRegion, nAge, length(a))) # P(Atilda = alow | W = w, A, T = 1)
Patildalow[,,1]     = rep(0.9^(seq(nAge) / 10), each = 3)
Patildalow[,,2]     = 0.01
Patildamedium       = array(0, c(nRegion, nAge, length(a))) # P(Atilda = amedium | W = w, A, T = 1)
Patildamedium[,,1]  = 0.2
Patildamedium[,,2]  = rep(0.9^(seq(nAge) / 10), each = 3)
Patildamedium[,,3]  = 0.01
Patildahigh         = array(0, c(nRegion, nAge, length(a))) # P(Atilda = ahigh | W = w, A, T = 1)
Patildahigh[,,1]    = 0.02
Patildahigh[,,2]    = 0.2
Patildahigh[,,3]    = rep(0.99^(seq(nAge) / 10), each = 3)

# P(Atilda | W, A, T)
Patilda         = array(0, c(wRegion = nRegion, wAge = nAge, a = length(a), atilda = length(a), t = 2))
Ptildatot       = Patildalow + Patildamedium + Patildahigh
Patilda[,,,1,2] = Patildalow / Ptildatot # P(Atilda | W, A, T = 1)
Patilda[,,,2,2] = Patildamedium / Ptildatot
Patilda[,,,3,2] = Patildahigh / Ptildatot
for(p in seq(length(a))) {
  Patilda[,,p,p,1] = 1 # we assume P(Atilda = alow | W = w, A = alow, T = 0) = 1
}

pY2       = c(0.07, 0.03, 0.0145, 0.055, 0.075, 0.04, 0.02, 0.055, 0.065, 
              0.075, 0.09, 0.03, 0.095, 0.02, 0.09, 0.07, 0.04, 0.07, 0.085,
              0.06) # incidence of the 20 non-targeted HPV types
NNT       = length(pY2) # 20 non-targeted HPV types

beta_1    = - 0.73  # direct (immunological) effect of the vaccine on Y1
beta_2    = 0       # the vaccine does not have an immunological effect on the 
# HPV viruses that are not targeted by the vaccine

# parameter values used for the simulation of Y_1 and Y2^(j) given T, A and W
q_WAge    = 0.01                # direct effect of WAge on the targeted type
q_WRegion = param$q_WRegion[,2] # direct effect of WRegion on the targeted type         
s_WAge    = param$s_WAge        # direct effect of WAge on each nontargeted type
s_WRegion = param$s_WRegion     # direct effect of WRegion on each nontargeted type

Nreplic   = 5*10^3 # number of study replications

### Function to be run for each scenario ---------------------------------------
Onerun = function(p){ 
  set.seed(1234)
  
  ## Parameter values in the scenario ------------------------------------------
  pY1       = PARAM[p, 1]
  a         = PARAM[p, 2:4]
  a_low     = a[1]
  a_medium  = a[2]
  a_high    = a[3]
  n         = PARAM[p, 5]
  
  ## Generation of the region of the n individuals -----------------------------
  WWRegion  = sample(wRegion, n, replace = TRUE, prob = PwRegion)        
  
  ## Generation of the age of the n individuals --------------------------------
  WWAge       = rep(NA, n)
  for(i in seq(nRegion)){
    id        = which(WWRegion == wRegion[i])
    WWAge[id] = sample(wAge, length(id), replace = TRUE, prob = PwAge[i,])
  }
  
  ## Computation of the dataset characteristics, analytically ------------------
  
  # Mean and variance of the unobserved confounder A
  meanA = 0
  varA  = 0
  for(i in seq(nRegion)){
    for(j in seq(nAge)){
      pw    = PwAge[i,j] * PwRegion[i]
      meanA = meanA + sum(a * Pa[i,j,]) * pw
      varA  = varA + sum(a^2 * Pa[i,j,]) * pw
    }}
  varA = varA - meanA^2

  # Mean of the treatment T, i.e., P(T = 1), and computation of P(T = 1 | W, A)
  pT    = array(NA, c(nRegion, nAge, length(a))) # P(T = 1 | W, A)
  meanT = 0
  for(i in seq(nRegion)){
    for(j in seq(nAge)){
      for(p in 1:length(a)){
        pt_aw     = plogis(alpha_T + delta_T * wAge[j] + gamma_T * wRegion[i] + zeta_T * a[p]) # logistic function
        pT[i,j,p] = pt_aw
        meanT     = meanT + pt_aw * Pa[i,j,p] * PwAge[i,j] * PwRegion[i]
      }}}
  
  # Mean and variance of the unobserved mediator Atilda
  meanAtilda  = 0
  varAtilda   = 0
  for(i in seq(nRegion)){
    for(j in seq(nAge)){
      for(p in 1:length(a)){
        pt_aw       = pT[i,j,p]
        meanAtilda  = meanAtilda + sum(a * (Patilda[i,j,p,,1] * (1 - pt_aw) + Patilda[i,j,p,,2] * pt_aw)) * Pa[i,j,p] * PwAge[i,j] * PwRegion[i] # sum t, a, w atilda P(Atilda | W, A, T) P(T | W, A) P(A | W) P(W) with P(T = O | W, A) = 1 - P(T = 1 | W, A)
        varAtilda   = varAtilda + sum(a^2 * (Patilda[i,j,p,,1] * (1 - pt_aw) + Patilda[i,j,p,,2] * pt_aw)) * Pa[i,j,p] * PwAge[i,j] * PwRegion[i]
      }}} 
  varAtilda         = varAtilda - meanAtilda^2
  
  # Mean of the targeted (i.e., primary) outcome Y1, and mean of the untargeted (individual) outcomes Y2
  meanY1 = 0
  meanY2 = 0
  for(i in seq(nRegion)){
    for(j in seq(nAge)){
      for(p in 1:length(a)){
        pt_aw   = pT[i,j,p]
        paw     = Pa[i,j,p] * PwAge[i,j] * PwRegion[i]
        exp_1   = exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i]))
        exp_2   = exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))
        meanY1  = meanY1 + exp_1 * sum(a * ((1 - pt_aw) * Patilda[i,j,p,,1] + pt_aw * Patilda[i,j,p,,2] * exp(beta_1))) * paw
        meanY2  = meanY2 + exp_2 * sum(a * ((1 - pt_aw) * Patilda[i,j,p,,1] + pt_aw * Patilda[i,j,p,,2])) * paw
      }}}

  # Computation of the intercepts, so that E(Y1) = pY1 and E(Y2) = pY2 for each HPV type.
  alpha_1         = log(pY1 / meanY1) # so that E(Y1) = pY1, for the targeted type
  alpha_2         = log(pY2 / meanY2) # so that E(Y2) = pY2, for each non-targeted type
  Alpha_1         = matrix(rep(alpha_1,n), nrow = n, byrow = TRUE)
  Alpha_2         = matrix(rep(alpha_2,n), nrow = n, byrow = TRUE)
  
  # Mean and variance of the targeted (i.e., primary) outcome Y1
  meanY1  = pY1
  varY1   = meanY1 * (1 - meanY1)
  
  # Mean of the untargeted (i.e., secondary) outcome Y2s
  meanY2s = sum(pY2) # E(Y2s)
  
  # Covariance between targeted and untargeted (i.e., primary and secondary) outcomes
  covY1Y2s  = 0 # cov(Y1c, Y2s)
  covY2     = 0 # cov(Y2^k,Y2^l)
  for(i in seq(nRegion)){
    for(j in seq(nAge)){
      for(p in 1:length(a)){
        pt_aw     = pT[i,j,p]
        paw       = Pa[i,j,p] * PwAge[i,j] * PwRegion[i]
        exp_1     = exp(alpha_1 + wAge[j] %*% t(q_WAge) + t(q_WRegion[i]))
        exp_2     = exp(alpha_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))
        covY1Y2s  = covY1Y2s + t(exp_2) %*% exp_1 * sum(((1 - pt_aw) * Patilda[i,j,p,,1] + pt_aw * Patilda[i,j,p,,2] * exp(beta_1)) * a^2) * paw
        covY2     = covY2 + t(exp_2) %*% exp_2 * sum(a^2 * ((1 - pt_aw) * Patilda[i,j,p,,1] + pt_aw * Patilda[i,j,p,,2]) * paw)
      }}}
  covY1Y2s = sum(covY1Y2s) - meanY1 * meanY2s 
  
  # Mean and variance of the untargeted (i.e., secondary) outcome Y2s
  varY2s  = meanY2s * (1 - meanY2s) + (sum(covY2 - diag(diag(covY2))))
  
  # Correlation between targeted and untargeted outcomes
  corrY1Y2s = covY1Y2s / sqrt(varY1 * varY2s)
  
  # Parameters of potential causal interest
  # log-ATE, log-NDE and log-NIE
  EY1_T1 = 0
  EY1_T0 = 0
  EY1_T1AtildaT0 = 0
  for(i in seq(nRegion)){
    for(j in seq(nAge)){
      for(p in 1:length(a)){
        paw     = Pa[i,j,p] * PwAge[i,j] * PwRegion[i]
        exp_1   = exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i]))
        EY1_T1  = EY1_T1 + exp_1 * sum(a * Patilda[i,j,p,,2]) * paw
        EY1_T0  = EY1_T0 + exp_1 * sum(a * Patilda[i,j,p,,1]) * paw
        EY1_T1AtildaT0  = EY1_T1AtildaT0 + exp_1 * sum(a * Patilda[i,j,p,,1]) * paw
      }}}
  EY1_T1          = EY1_T1 * exp(alpha_1 + beta_1)
  EY1_T0          = EY1_T0 * exp(alpha_1)
  EY1_T1AtildaT0  = EY1_T1AtildaT0 * exp(alpha_1 + beta_1)
  
  log.ATE.true  = log(EY1_T1 / EY1_T0) # overall log-ratio effect of the vaccine
  log.NIE.true  = log(EY1_T1 / EY1_T1AtildaT0) # harmful indirect (behavioral) log-ratio effect of the vaccine
  log.NDE.true  = log(EY1_T1AtildaT0 / EY1_T0) # protective direct (immunological) log-ratio effect of the vaccine
  
  res           = NULL
  
  for(nrep in 1:Nreplic){
    
    ## Data generation ---------------------------------------------------------
    sample  = sample(n, n)
    
    WRegion = WWRegion[sample]
    WAge    = WWAge[sample]
    
    A = rep(NA, n)
    for(i in seq(nRegion)){
      for(j in seq(nAge)){
        id    = which((WRegion == wRegion[i])&(WAge == wAge[j]))
        A[id] = sample(c(a_low, a_medium, a_high), size = length(id), 
                       replace = TRUE, prob = Pa[i,j,])
      }}
    
    T = rbinom(n, size = 1, prob = exp(alpha_T + delta_T * WAge + gamma_T * WRegion + zeta_T * A) / (1 + exp(alpha_T + delta_T * WAge + gamma_T * WRegion + zeta_T * A)))
    
    Atilda = rep(NA, n)
    Atilda[T == 0] = A[T == 0]
    for(i in seq(nRegion)){
      for(j in seq(nAge)){
        for(p in 1:length(a)){
          id          = which((WRegion == wRegion[i])&(WAge == wAge[j])&(T == 1)&(A == a[p]))
          Atilda[id]  = sample(c(a_low, a_medium, a_high), size = length(id), 
                               replace = TRUE, prob = Patilda[i,j,p,,2])
      }}}
    
    WRegion_Y1 = 0
    for(i in seq(nRegion)){
      WRegion_Y1 = WRegion_Y1 + 1 * (WRegion == wRegion[i]) %*% t(q_WRegion[i])
    }
    
    p1              = exp(Alpha_1 + T %*% t(beta_1) + WAge %*% t(q_WAge) + WRegion_Y1)
    PY1             = Atilda * p1
    Y1              = rbinom(n, size = 1, prob = PY1) 
    
    p2              = exp(Alpha_2 + WAge %*% t(s_WAge) + WRegion %*% t(s_WRegion))
    PY2             = Atilda * p2 
    Y2              = matrix(rbinom(NNT*n, size = 1, prob = PY2), nrow = n, byrow = FALSE) 
    Y2s             = rowSums(Y2)   

    W = interaction(WAge, WRegion, sep = "_")
    
    ## Estimation with the Joint methods ---------------------------------------
    
    est.jointReg    = JointReg(Y1 = Y1, Y2 = Y2s, T = T, Wcont = WAge, 
                               Wcat = WRegion) # estimation of the treatment effect on Y1 using Joint-Reg
    est.jointMH     = JointMH(Y1 = Y1, Y2 = Y2s, T = T, W = W) # estimation of the treatment effect on Y1 using Joint-MH
    est.jointNC     = JointNC(Y1 = Y1, Y2 = Y2s, T = T) # estimation of the treatment effect on Y1 using Joint-NC (i.e., when the observed confounders are not used)
 
    jointReg        = cbind("Joint-Reg", est.jointReg$beta_1.hat, 
                            est.jointReg$se.beta_1.hat, 
                            est.jointReg$beta_1.hat - est.jointReg$se.beta_1.hat * qnorm(0.975), 
                            est.jointReg$beta_1.hat + est.jointReg$se.beta_1.hat * qnorm(0.975))
    jointMH         = cbind("Joint-MH", est.jointMH$beta_1.hat, 
                            est.jointMH$se.beta_1.hat,  
                            est.jointMH$beta_1.hat - est.jointMH$se.beta_1.hat * qnorm(0.975), 
                            est.jointMH$beta_1.hat + est.jointMH$se.beta_1.hat * qnorm(0.975))
    jointNC         = cbind("Joint-NC", est.jointNC$beta_1.hat, 
                            est.jointNC$se.beta_1.hat, 
                            est.jointNC$beta_1.hat - est.jointNC$se.beta_1.hat * qnorm(0.975), 
                            est.jointNC$beta_1.hat + est.jointNC$se.beta_1.hat * qnorm(0.975))
    MH              = cbind("MH", est.jointMH$beta_1.hat.naive, 
                            est.jointMH$se.beta_1.hat.naive, 
                            est.jointMH$beta_1.hat.naive - est.jointMH$se.beta_1.hat.naive * qnorm(0.975), 
                            est.jointMH$beta_1.hat.naive + est.jointMH$se.beta_1.hat.naive * qnorm(0.975)) # naive MH estimation of the treatment effect on Y1 using, without using Y2
    Reg             = cbind("Reg", est.jointReg$beta_1.hat.naive, 
                            est.jointReg$se.beta_1.hat.naive, 
                            est.jointReg$beta_1.hat.naive - est.jointReg$se.beta_1.hat.naive * qnorm(0.975), 
                            est.jointReg$beta_1.hat.naive + est.jointReg$se.beta_1.hat.naive * qnorm(0.975)) # naive Reg estimation of the treatment effect on Y1 using, without using Y2
 
    recap = rbind(jointReg, jointMH, jointNC, MH, Reg)
    colnames(recap) = c("Approach", "beta_1.hat", "se.beta_1.hat", "CI.left", "CI.right")
    cbind(recap) 

    res = rbind(res, cbind(recap, n = n, pY1 = pY1, 
                           pY2 = paste(pY2, collapse = "_"), 
                           A = paste(a, collapse = "_"), 
                           log.ATE.true = as.numeric(log.ATE.true), 
                           log.NDE.true = as.numeric(log.NDE.true), 
                           log.NIE.true = as.numeric(log.NIE.true), 
                           corrY1Y2s.true = as.numeric(corrY1Y2s), 
                           corrY1Y2s = cor(Y1,Y2s), 
                           meanY1.true = as.numeric(meanY1), meanY1 = mean(Y1), 
                           meanY2s.true = meanY2s, meanY2s = mean(Y2s), 
                           meanA.true = meanA, meanA = mean(A), 
                           varA.true = varA, varA = var(A), 
                           meanAtilda.true = meanAtilda, 
                           meanAtilda = mean(Atilda), 
                           varAtilda.true = varAtilda, varAtilda = var(Atilda)))
  }
  myfile  = paste0("RES_Observational-n", n, "-pY1", pY1, "-A", 
                   paste(a, collapse = "_"), "-beta1", 
                   round(beta_1, digits = 3), ".RData")
  save(res, file = myfile)
}
 
mclapply(1:nrow(PARAM), Onerun, mc.cores = 18)

### Combining the results from the different scenarios -------------------------
RES = NULL
for(p in 1: nrow(PARAM)){
  pY1             = PARAM[p, 1]
  a               = PARAM[p, 2:4]
  n               = PARAM[p, 5]
  myfile          = paste0("RES_Observational-n", n, "-pY1", pY1, "-A", 
                           paste(a, collapse = "_"), "-beta1", 
                           round(beta_1, digits = 3), ".RData")
  load(myfile)
  RES = rbind(RES, res)
}
RECAP           = as.data.frame(RES)
ColNames        = colnames(RECAP[,c(2:7, 10:26)])
RECAP[ColNames] = sapply(RECAP[ColNames], as.numeric)
RECAP$pY1       = as.factor(RECAP$pY1)
RECAP$n         = as.factor(RECAP$n)
RECAP$A         = as.factor(RECAP$A)
RECAP$Approach  = as.factor(RECAP$Approach)

myfile  = paste0("RECAP_Observational-beta1", round(beta_1, digits = 3), ".RData")
save(RECAP, file = myfile)

## Analysis of the simulation results ------------------------------------------
Res = NULL
for(i in 1:nrow(PARAM)){
  RECAP1 = RECAP[((i-1)*(5*Nreplic) + 1):(i*(5*Nreplic)),]

  # Coverage of the 95% confidence intervals
  cov.JointNC     = sum((RECAP1[which(RECAP1$Approach == "Joint-NC"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Joint-NC"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Joint-NC"),5])
  cov.JointMH     = sum((RECAP1[which(RECAP1$Approach == "Joint-MH"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Joint-MH"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Joint-MH"),5])
  cov.MH          = sum((RECAP1[which(RECAP1$Approach == "MH"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "MH"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "MH"),5])
  cov.JointReg    = sum((RECAP1[which(RECAP1$Approach == "Joint-Reg"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Joint-Reg"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Joint-Reg"),5])
  cov.Reg         = sum((RECAP1[which(RECAP1$Approach == "Reg"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Reg"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Reg"),5])

  # Standard deviation over the 5000 replications
  sd.JointNC      = sd(RECAP1[which(RECAP1$Approach == "Joint-NC"),2])
  sd.JointMH      = sd(RECAP1[which(RECAP1$Approach == "Joint-MH"),2])
  sd.MH           = sd(RECAP1[which(RECAP1$Approach == "MH"),2])
  sd.JointReg     = sd(RECAP1[which(RECAP1$Approach == "Joint-Reg"),2])
  sd.Reg          = sd(RECAP1[which(RECAP1$Approach == "Reg"),2])

  # Mean sandwich standard error over the 5000 replications
  sandwich_se.JointNC   = mean(RECAP1[which(RECAP1$Approach == "Joint-NC"),3])
  sandwich_se.JointMH   = mean(RECAP1[which(RECAP1$Approach == "Joint-MH"),3])
  sandwich_se.MH        = mean(RECAP1[which(RECAP1$Approach == "MH"),3])
  sandwich_se.JointReg  = mean(RECAP1[which(RECAP1$Approach == "Joint-Reg"),3])
  sandwich_se.Reg       = mean(RECAP1[which(RECAP1$Approach == "Reg"),3])

  # Mean effect estimate and mean squared error over the 5000 replications
  mean.JointNC    = mean(RECAP1[which(RECAP1$Approach == "Joint-NC"),2])
  MSE.JointNC     = mean((RECAP1[which(RECAP1$Approach == "Joint-NC"),2] - beta_1)^2)
  mean.JointMH    = mean(RECAP1[which(RECAP1$Approach == "Joint-MH"),2])
  MSE.JointMH     = mean((RECAP1[which(RECAP1$Approach == "Joint-MH"),2] - beta_1)^2)
  mean.MH         = mean(RECAP1[which(RECAP1$Approach == "MH"),2])
  MSE.MH          = mean((RECAP1[which(RECAP1$Approach == "MH"),2] - beta_1)^2)
  mean.JointReg   = mean(RECAP1[which(RECAP1$Approach == "Joint-Reg"),2])
  MSE.JointReg    = mean((RECAP1[which(RECAP1$Approach == "Joint-Reg"),2] - beta_1)^2)
  mean.Reg        = mean(RECAP1[which(RECAP1$Approach == "Reg"),2])
  MSE.Reg         = mean((RECAP1[which(RECAP1$Approach == "Reg"),2] - beta_1)^2)

  # Ratio of MSE for the naive method MH and Joint methods
  eff_MH_JointMH    = MSE.MH / MSE.JointMH
  eff_MH_JointNC    = MSE.MH / MSE.JointNC
  eff_MH_JointReg   = MSE.MH / MSE.JointReg

  # Ratio of MSE for the naive method Reg and Joint methods
  eff_Reg_JointMH   = MSE.Reg / MSE.JointMH
  eff_Reg_JointNC   = MSE.Reg / MSE.JointNC
  eff_Reg_JointReg  = MSE.Reg / MSE.JointReg

  Res = rbind(Res, c(eff_MH_JointMH = eff_MH_JointMH,
                     eff_MH_JointReg = eff_MH_JointReg,
                     eff_MH_JointNC = eff_MH_JointNC,
                     eff_Reg_JointMH = eff_Reg_JointMH,
                     eff_Reg_JointReg = eff_Reg_JointReg,
                     eff_Reg_JointNC = eff_Reg_JointNC,
                     relbias.JointMH = abs((mean.JointMH - beta_1) / beta_1),
                     relbias.JointReg = abs((mean.JointReg - beta_1) / beta_1),
                     relbias.JointNC = abs((mean.JointNC - beta_1) / beta_1),
                     relbias.MH = abs((mean.MH - beta_1) / beta_1),
                     relbias.Reg = abs((mean.Reg - beta_1) / beta_1),
                     empir_sd.JointMH = sd.JointMH,
                     sandwich_se.JointMH = sandwich_se.JointMH,
                     empir_sd.JointReg = sd.JointReg,
                     sandwich_se.JointReg = sandwich_se.JointReg,
                     empir_sd.JointNC = sd.JointNC,
                     sandwich_se.JointNC = sandwich_se.JointNC,
                     empir_sd.MH = sd.MH,
                     sandwich_se.MH = sandwich_se.MH,
                     empir_sd.Reg = sd.Reg,
                     sandwich_se.Reg = sandwich_se.Reg,
                     CIcov.JointMH = cov.JointMH, CIcov.JointReg = cov.JointReg,
                     CIcov.JointNC = cov.JointNC,
                     CIcov.MH = cov.MH, CIcov.Reg = cov.Reg,
                     n = as.character(RECAP1[1,]$n),
                     pY1 = as.character(RECAP1[1,]$pY1),
                     A =  as.character(RECAP1[1,]$A),
                     beta_1 = beta_1,
                     log.ATE.true = RECAP1[1,]$log.ATE.true,
                     log.NDE.true = RECAP1[1,]$log.NDE.true,
                     log.NIE.true = RECAP1[1,]$log.NIE.true,
                     corrY1Y2s.true = RECAP1[1,]$corrY1Y2s.true,
                     corrY1Y2s = mean(RECAP1$corrY1Y2s),
                     meanY1.true = RECAP1[1,]$meanY1.true,
                     meanY1 = mean(RECAP1$meanY1),
                     meanY2s.true = RECAP1[1,]$meanY2s.true,
                     meanY2s = mean(RECAP1$meanY2s),
                     meanA.true = RECAP1[1,]$meanA.true,
                     meanA = mean(RECAP1$meanA),
                     varA.true = RECAP1[1,]$varA.true,
                     varA = mean(RECAP1$varA),
                     meanAtilda.true = RECAP1[1,]$meanAtilda.true,
                     meanAtilda = mean(RECAP1$meanAtilda),
                     varAtilda.true = RECAP1[1,]$varAtilda.true,
                     varAtilda = mean(RECAP1$varAtilda)
  ))
}
Res           = as.data.frame(Res)
ColNames      = colnames(Res[,c(1:27, 30:47)])
Res[ColNames] = sapply(Res[ColNames], as.numeric)
save(Res, file = paste0("Res_Observational-beta1", round(beta_1, digits = 3),
                        ".RData"))
Res[ColNames] = round(Res[ColNames], digits = 3)
write.csv(Res, file = paste0("Res_Observational-beta1",
                             round(beta_1, digits = 3), ".csv"))
write.xlsx(Res, file = paste0("Res_Observational-beta1",
                              round(beta_1, digits = 3), ".xlsx"))
Res1 = Res[which(Res$n == 10000),]
latextable = cbind(Res1[, c(7,8,10,11,12,14,13,15,28,29,34)]) # only the scenarios with n = 10,000
print(xtable(latextable, digits = 3), include.rownames=FALSE) # print the latex table
Res2 = Res[which(Res$n == 5000),]
latextable2 = cbind(Res2[, c(7,8,10,11,12,14,13,15,28,29,34)])
print(xtable(latextable2, digits = 3), include.rownames=FALSE) # latex table for n = 5,000

## Plotting the results --------------------------------------------------------
RECAP$pY1 = factor(RECAP$pY1,
                   labels = c(expression(P(Y[1] == 1) == 0.025),
                              expression(P(Y[1] == 1) == 0.05),
                              expression(P(Y[1] == 1) == 0.14)))
RECAP$A = factor(RECAP$A,
                 labels = c(expression(a[low]==~0~","~a[medium]==~0.75~","~a[high]==~1.5),
                            expression(a[low]==~0~","~a[medium]==~1~","~a[high]==~2),
                            expression(a[low]==~0~","~a[medium]==~1~","~a[high]==~2.5) ))
RECAP$Approach = factor(RECAP$Approach,levels(RECAP$Approach)[c(1,3,2,4,5)])
RECAP = RECAP[-which(RECAP$Approach == "Joint-NC"),]
plot = ggplot(RECAP, aes(x = n, y = beta_1.hat, color = Approach)) +
  geom_boxplot(coef = NULL) + geom_hline(aes(yintercept = log.NDE.true)) +
  geom_hline(aes(yintercept = log.ATE.true), linetype = "dotdash") + theme_light() +
  ylab("Estimated log-ratio effect of the vaccine on the targeted outcome") +
  xlab("Observational study size") +
  facet_grid(pY1~A, labeller = label_parsed, scales = "free_y") +
  theme(plot.title = element_text(size = 11),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 11),
        strip.background = element_rect(color="black", fill="white", size = 0.5, linetype="solid"),
        strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black"))
pdf(paste0("Comparison_Observational-beta1", round(beta_1, digits = 3), ".pdf"),
    width = 10, height = 7)
plot
dev.off()

