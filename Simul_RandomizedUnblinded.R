## -----------------------------------------------------------------------------
## Description: This script illustrates estimation of the NIE when using the 
##              methods proposed in Etievant et al. (Biometrics, 2023) in an 
##              randomized setting with risk compensation
##
##              The simulation is displayed in Web Appendix ?? in the 
##              Supplementary Material of Etievant (2025)
##
##              Comments:
##              Treatment (i.e., vaccine) assignment has been randomized
##              The vaccine has a direct and an indirect effect on the targeted 
##              outcome because of risk compensation. The mediator is unmeasured
## -----------------------------------------------------------------------------

### load packages and files ----------------------------------------------------
source("EstimationFunctions.R")
load("Parameters.RData") # file with certain of the fixed parameters
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
#   A       unobserved covariate e.g. sexual activity pre-vaccination
#
#   W       observed covariate e.g. age and geographical region
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

pT = 0.5 # probability of receiving the vaccine is 0.5

a                   = Aa[1, ] # Atilda can also take values alow, amedium and 
# ahigh, with probability depending on the values of W, A and T
Patildalow          = array(0, c(nRegion, nAge, length(a))) # P(Atilda = alow | W = w, A, T = 1)
Patildalow[,,1]     = rep(0.9^(seq(nAge) / 10), each = 3)
Patildalow[,,2]     = 0.01
Patildamedium       = array(0, c(nRegion, nAge, length(a))) # P(Atilda = amedium | W = w, A, T = 1)
Patildamedium[,,1]  = 0.08
Patildamedium[,,2]  = rep(0.9^(seq(nAge) / 10), each = 3)
Patildamedium[,,3]  = 0.01
Patildahigh         = array(0, c(nRegion, nAge, length(a))) # P(Atilda = ahigh | W = w, A, T = 1)
Patildahigh[,,1]    = 0.02
Patildahigh[,,2]    = 0.1
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
  
  # Mean and variance of the unobserved mediator Atilda
  meanAtilda  = 0
  varAtilda   = 0
  for(i in seq(nRegion)){
    for(j in seq(nAge)){
      for(p in 1:length(a)){
        meanAtilda  = meanAtilda + sum(a * (Patilda[i,j,p,,1] * (1 - pT) + Patilda[i,j,p,,2] * pT)) * Pa[i,j,p] * PwAge[i,j] * PwRegion[i] # sum t, a, w atilda P(Atilda | W, A, T) P(T) P(A | W) P(W) with P(T = O | W, A) = 1 - P(T = 1 | W, A)
        varAtilda   = varAtilda + sum(a^2 * (Patilda[i,j,p,,1] * (1 - pT) + Patilda[i,j,p,,2] * pT)) * Pa[i,j,p] * PwAge[i,j] * PwRegion[i]
      }}} 
  varAtilda         = varAtilda - meanAtilda^2
  
  # Mean of the targeted (i.e., primary) outcome Y1, and mean of the untargeted (individual) outcomes Y2
  meanY1 = 0
  meanY2 = 0
  for(i in seq(nRegion)){
    for(j in seq(nAge)){
      for(p in 1:length(a)){
        paw     = Pa[i,j,p] * PwAge[i,j] * PwRegion[i]
        exp_1   = exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i]))
        exp_2   = exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))
        meanY1  = meanY1 + exp_1 * sum(a * ((1 - pT) * Patilda[i,j,p,,1] + pT * Patilda[i,j,p,,2] * exp(beta_1))) * paw
        meanY2  = meanY2 + exp_2 * sum(a * ((1 - pT) * Patilda[i,j,p,,1] + pT * Patilda[i,j,p,,2])) * paw
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
        paw       = Pa[i,j,p] * PwAge[i,j] * PwRegion[i]
        exp_1     = exp(alpha_1 + wAge[j] %*% t(q_WAge) + t(q_WRegion[i]))
        exp_2     = exp(alpha_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))
        covY1Y2s  = covY1Y2s + t(exp_2) %*% exp_1 * sum(((1 - pT) * Patilda[i,j,p,,1] + pT * Patilda[i,j,p,,2] * exp(beta_1)) * a^2) * paw
        covY2     = covY2 + t(exp_2) %*% exp_2 * sum(a^2 * ((1 - pT) * Patilda[i,j,p,,1] + pT * Patilda[i,j,p,,2]) * paw)
      }}}
  covY1Y2s = sum(covY1Y2s) - meanY1 * meanY2s 
  
  # Mean and variance of the untargeted (i.e., secondary) outcome Y2s
  varY2s  = meanY2s * (1 - meanY2s) + (sum(covY2 - diag(diag(covY2))))
  
  # Correlation between targeted and untargeted outcomes
  corrY1Y2s = covY1Y2s / sqrt(varY1 * varY2s)
 
  # Parameters of potential causal interest
  # ATE, NDE and NIE
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
  
  ATE.true  = log(EY1_T1 / EY1_T0) # overall protective of the vaccine
  NIE.true  = log(EY1_T1 / EY1_T1AtildaT0) # harmful indirect (behavioral) effect of the vaccine
  NDE.true  = log(EY1_T1AtildaT0 / EY1_T0) # protective direct (immunological) effect of the vaccine
  
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
    
    T = rbinom(n, size = 1, prob = 1/2)
    
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
                               Wcat = WRegion) # estimation with Joint-Reg
    est.jointMH     = JointMH(Y1, Y2s, T = T, W = W) # estimation with method Joint-MH
    est.jointNC     = JointNC(Y1 = Y1, Y2 = Y2s, T = T) # estimation with method Joint-NC (i.e., when the observed covariates are not used)
 
    Reg.NDE     = cbind("Reg.NDE", est.jointReg$beta_1.hat, 
                        est.jointReg$se.beta_1.hat, 
                        est.jointReg$beta_1.hat - est.jointReg$se.beta_1.hat * qnorm(0.975), 
                        est.jointReg$beta_1.hat + est.jointReg$se.beta_1.hat * qnorm(0.975)) # beta1 - beta2 should gives- NDE
    Reg.ATE     = cbind("Reg.ATE", est.jointReg$beta_1.hat.naive, 
                        est.jointReg$se.beta_1.hat.naive, 
                        est.jointReg$beta_1.hat.naive - est.jointReg$se.beta_1.hat.naive * qnorm(0.975), 
                        est.jointReg$beta_1.hat.naive + est.jointReg$se.beta_1.hat.naive * qnorm(0.975)) # beta1 should give ATE
    Reg.NIE     = cbind("Reg.NIE", est.jointReg$beta_2.hat.naive, 
                        est.jointReg$se.beta_2.hat.naive, 
                        est.jointReg$beta_2.hat.naive - est.jointReg$se.beta_2.hat.naive * qnorm(0.975), 
                        est.jointReg$beta_2.hat.naive + est.jointReg$se.beta_2.hat.naive * qnorm(0.975)) # beta2 should give NIE
    
    MH.NDE      = cbind("MH.NDE", est.jointMH$beta_1.hat, 
                        est.jointMH$se.beta_1.hat,  
                        est.jointMH$beta_1.hat - est.jointMH$se.beta_1.hat * qnorm(0.975), 
                        est.jointMH$beta_1.hat + est.jointMH$se.beta_1.hat * qnorm(0.975))
    MH.ATE      = cbind("MH.ATE", est.jointMH$beta_1.hat.naive, 
                        est.jointMH$se.beta_1.hat.naive, 
                        est.jointMH$beta_1.hat.naive - est.jointMH$se.beta_1.hat.naive * qnorm(0.975), 
                        est.jointMH$beta_1.hat.naive + est.jointMH$se.beta_1.hat.naive * qnorm(0.975))
    MH.NIE      = cbind("MH.NIE", est.jointMH$beta_2.hat.naive, 
                        est.jointMH$se.beta_2.hat.naive, 
                        est.jointMH$beta_2.hat.naive - est.jointMH$se.beta_2.hat.naive * qnorm(0.975), 
                        est.jointMH$beta_2.hat.naive + est.jointMH$se.beta_2.hat.naive * qnorm(0.975))
    
    NC.NDE      = cbind("NC.NDE", est.jointNC$beta_1.hat, 
                        est.jointNC$se.beta_1.hat, 
                        est.jointNC$beta_1.hat - est.jointNC$se.beta_1.hat * qnorm(0.975), 
                        est.jointNC$beta_1.hat + est.jointNC$se.beta_1.hat * qnorm(0.975))
    NC.ATE      = cbind("NC.ATE",  est.jointNC$beta_1.hat.naive, 
                        est.jointNC$se.beta_1.hat.naive, 
                        est.jointNC$beta_1.hat.naive - est.jointNC$se.beta_1.hat.naive * qnorm(0.975), 
                        est.jointNC$beta_1.hat.naive + est.jointNC$se.beta_1.hat.naive * qnorm(0.975))
    NC.NIE      = cbind("NC.NIE", est.jointNC$beta_2.hat.naive, 
                        est.jointNC$se.beta_2.hat.naive, 
                        est.jointNC$beta_2.hat.naive - est.jointNC$se.beta_2.hat.naive * qnorm(0.975), 
                        est.jointNC$beta_2.hat.naive + est.jointNC$se.beta_2.hat.naive * qnorm(0.975))
    
    recap = rbind(Reg.NDE, Reg.NIE, Reg.ATE, MH.NDE, MH.NIE, MH.ATE, NC.NDE, 
                  NC.NIE, NC.ATE)
    colnames(recap) = c("Approach", "Effect.hat", "se.Effect.hat", 
                        "CI.Effect.left", "CI.Effect.right")
    
    res = rbind(res, cbind(recap, n = n, pY1 = pY1, 
                           pY2 = paste(pY2, collapse = "_"), 
                           A = paste(a, collapse = "_"), 
                           ATE.true = as.numeric(ATE.true), 
                           NDE.true = as.numeric(NDE.true), 
                           NIE.true = as.numeric(NIE.true), 
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
  myfile  = paste0("RES_RandomizedUnblinded-n", n, "-pY1", pY1, "-A", 
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
  myfile  = paste0("RES_RandomizedUnblinded-n", n, "-pY1", pY1, "-A", 
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
myfile          = paste0("RECAP_RandomizedUnblinded-beta1", 
                         round(beta_1, digits = 3), ".RData")
save(RECAP, file = myfile)

## Analysis of the simulation results ------------------------------------------
# Details of the results
Nreplic   = 5*10^3
Res = NULL
for(i in 1:nrow(PARAM)){
  RECAP1 = RECAP[((i-1)*(9*Nreplic) + 1):(i*(9*Nreplic)),]

  ATE = RECAP1$ATE.true[1]
  NDE = RECAP1$NDE.true[1]
  NIE = RECAP1$NIE.true[1]

  # Coverage of the confidence intervals
  cov.ATE.NC     = sum((RECAP1[which(RECAP1$Approach == "NC.ATE"),4] < ATE)&(RECAP1[which(RECAP1$Approach == "NC.ATE"),5] > ATE)) / length(RECAP1[which(RECAP1$Approach == "NC.ATE"),5])
  cov.NDE.NC     = sum((RECAP1[which(RECAP1$Approach == "NC.NDE"),4] < NDE)&(RECAP1[which(RECAP1$Approach == "NC.NDE"),5] > NDE)) / length(RECAP1[which(RECAP1$Approach == "NC.NDE"),5])
  cov.NIE.NC     = sum((RECAP1[which(RECAP1$Approach == "NC.NIE"),4] < NIE)&(RECAP1[which(RECAP1$Approach == "NC.NIE"),5] > NIE)) / length(RECAP1[which(RECAP1$Approach == "NC.NIE"),5])
  cov.ATE.MH     = sum((RECAP1[which(RECAP1$Approach == "MH.ATE"),4] < ATE)&(RECAP1[which(RECAP1$Approach == "MH.ATE"),5] > ATE)) / length(RECAP1[which(RECAP1$Approach == "MH.ATE"),5])
  cov.NDE.MH     = sum((RECAP1[which(RECAP1$Approach == "MH.NDE"),4] < NDE)&(RECAP1[which(RECAP1$Approach == "MH.NDE"),5] > NDE)) / length(RECAP1[which(RECAP1$Approach == "MH.NDE"),5])
  cov.NIE.MH     = sum((RECAP1[which(RECAP1$Approach == "MH.NIE"),4] < NIE)&(RECAP1[which(RECAP1$Approach == "MH.NIE"),5] > NIE)) / length(RECAP1[which(RECAP1$Approach == "MH.NIE"),5])
  cov.ATE.Reg     = sum((RECAP1[which(RECAP1$Approach == "Reg.ATE"),4] < ATE)&(RECAP1[which(RECAP1$Approach == "Reg.ATE"),5] > ATE)) / length(RECAP1[which(RECAP1$Approach == "Reg.ATE"),5])
  cov.NDE.Reg     = sum((RECAP1[which(RECAP1$Approach == "Reg.NDE"),4] < NDE)&(RECAP1[which(RECAP1$Approach == "Reg.NDE"),5] > NDE)) / length(RECAP1[which(RECAP1$Approach == "Reg.NDE"),5])
  cov.NIE.Reg     = sum((RECAP1[which(RECAP1$Approach == "Reg.NIE"),4] < NIE)&(RECAP1[which(RECAP1$Approach == "Reg.NIE"),5] > NIE)) / length(RECAP1[which(RECAP1$Approach == "Reg.NIE"),5])

  # Standard deviation over the 5000 replications
  sd.ATE.MH           = sd(RECAP1[which(RECAP1$Approach == "MH.ATE"),2])
  sd.NDE.MH           = sd(RECAP1[which(RECAP1$Approach == "MH.NDE"),2])
  sd.NIE.MH           = sd(RECAP1[which(RECAP1$Approach == "MH.NIE"),2])
  sd.ATE.NC           = sd(RECAP1[which(RECAP1$Approach == "NC.ATE"),2])
  sd.NDE.NC           = sd(RECAP1[which(RECAP1$Approach == "NC.NDE"),2])
  sd.NIE.NC           = sd(RECAP1[which(RECAP1$Approach == "NC.NIE"),2])
  sd.ATE.Reg           = sd(RECAP1[which(RECAP1$Approach == "Reg.ATE"),2])
  sd.NDE.Reg           = sd(RECAP1[which(RECAP1$Approach == "Reg.NDE"),2])
  sd.NIE.Reg           = sd(RECAP1[which(RECAP1$Approach == "Reg.NIE"),2])

  # Mean sandwich standard error over the 5000 replications
  sandwich_se.ATE.NC        = mean(RECAP1[which(RECAP1$Approach == "NC.ATE"),3])
  sandwich_se.NDE.NC        = mean(RECAP1[which(RECAP1$Approach == "NC.NDE"),3])
  sandwich_se.NIE.NC        = mean(RECAP1[which(RECAP1$Approach == "NC.NIE"),3])
  sandwich_se.ATE.MH        = mean(RECAP1[which(RECAP1$Approach == "MH.ATE"),3])
  sandwich_se.NDE.MH        = mean(RECAP1[which(RECAP1$Approach == "MH.NDE"),3])
  sandwich_se.NIE.MH        = mean(RECAP1[which(RECAP1$Approach == "MH.NIE"),3])
  sandwich_se.ATE.Reg        = mean(RECAP1[which(RECAP1$Approach == "Reg.ATE"),3])
  sandwich_se.NDE.Reg        = mean(RECAP1[which(RECAP1$Approach == "Reg.NDE"),3])
  sandwich_se.NIE.Reg        = mean(RECAP1[which(RECAP1$Approach == "Reg.NIE"),3])

  # Mean effect estimate and mean squared error over the 5000 replications
  mean.ATE.NC         = mean(RECAP1[which(RECAP1$Approach == "NC.ATE"),2])
  MSE.ATE.NC          = mean((RECAP1[which(RECAP1$Approach == "NC.ATE"),2] - ATE)^2) # MSE
  mean.NDE.NC         = mean(RECAP1[which(RECAP1$Approach == "NC.NDE"),2])
  MSE.NDE.NC          = mean((RECAP1[which(RECAP1$Approach == "NC.NDE"),2] - NDE)^2) # MSE
  mean.NIE.NC         = mean(RECAP1[which(RECAP1$Approach == "NC.NIE"),2])
  MSE.NIE.NC          = mean((RECAP1[which(RECAP1$Approach == "NC.NIE"),2] - NIE)^2) # MSE
  mean.ATE.MH         = mean(RECAP1[which(RECAP1$Approach == "MH.ATE"),2])
  MSE.ATE.MH          = mean((RECAP1[which(RECAP1$Approach == "MH.ATE"),2] - ATE)^2) # MSE
  mean.NDE.MH         = mean(RECAP1[which(RECAP1$Approach == "MH.NDE"),2])
  MSE.NDE.MH          = mean((RECAP1[which(RECAP1$Approach == "MH.NDE"),2] - NDE)^2) # MSE
  mean.NIE.MH         = mean(RECAP1[which(RECAP1$Approach == "MH.NIE"),2])
  MSE.NIE.MH          = mean((RECAP1[which(RECAP1$Approach == "MH.NIE"),2] - NIE)^2) # MSE
  mean.ATE.Reg         = mean(RECAP1[which(RECAP1$Approach == "Reg.ATE"),2])
  MSE.ATE.Reg          = mean((RECAP1[which(RECAP1$Approach == "Reg.ATE"),2] - ATE)^2) # MSE
  mean.NDE.Reg         = mean(RECAP1[which(RECAP1$Approach == "Reg.NDE"),2])
  MSE.NDE.Reg          = mean((RECAP1[which(RECAP1$Approach == "Reg.NDE"),2] - NDE)^2) # MSE
  mean.NIE.Reg         = mean(RECAP1[which(RECAP1$Approach == "Reg.NIE"),2])
  MSE.NIE.Reg          = mean((RECAP1[which(RECAP1$Approach == "Reg.NIE"),2] - NIE)^2) # MSE

  Res = rbind(Res, c(relbias.NDE.MH = abs((mean.NDE.MH - NDE) / NDE),
                     relbias.NIE.MH = abs((mean.NIE.MH - NIE) / NIE),
                     relbias.ATE.MH = abs((mean.ATE.MH - ATE) / ATE),
                     relbias.NDE.Reg = abs((mean.NDE.Reg - NDE) / NDE),
                     relbias.NIE.Reg = abs((mean.NIE.Reg - NIE) / NIE),
                     relbias.ATE.Reg = abs((mean.ATE.Reg - ATE) / ATE),
                     relbias.NDE.NC = abs((mean.NDE.NC - NDE) / NDE),
                     relbias.NIE.NC = abs((mean.NIE.NC - NIE) / NIE),
                     relbias.ATE.NC = abs((mean.ATE.NC - ATE) / ATE),
                     empir_sd.NDE.MH = sd.NDE.MH,
                     sandwich_se.NDE.MH = sandwich_se.NDE.MH,
                     empir_sd.NIE.MH = sd.NIE.MH,
                     sandwich_se.NIE.MH = sandwich_se.NIE.MH,
                     empir_sd.ATE.MH = sd.ATE.MH,
                     sandwich_se.ATE.MH = sandwich_se.ATE.MH,
                     empir_sd.NDE.Reg = sd.NDE.Reg,
                     sandwich_se.NDE.Reg = sandwich_se.NDE.Reg,
                     empir_sd.NIE.Reg = sd.NIE.Reg,
                     sandwich_se.NIE.Reg = sandwich_se.NIE.Reg,
                     empir_sd.ATE.Reg = sd.ATE.Reg,
                     sandwich_se.ATE.Reg = sandwich_se.ATE.Reg,
                     empir_sd.NDE.NC = sd.NDE.NC,
                     sandwich_se.NDE.NC = sandwich_se.NDE.NC,
                     empir_sd.NIE.NC = sd.NIE.NC,
                     sandwich_se.NIE.NC = sandwich_se.NIE.NC,
                     empir_sd.ATE.NC = sd.ATE.NC,
                     sandwich_se.ATE.NC = sandwich_se.ATE.NC,
                     CIcov.NDE.MH = cov.NDE.MH,
                     CIcov.NIE.MH = cov.NIE.MH,
                     CIcov.ATE.MH = cov.ATE.MH,
                     CIcov.NDE.Reg = cov.NDE.Reg,
                     CIcov.NIE.Reg = cov.NIE.Reg,
                     CIcov.ATE.Reg = cov.ATE.Reg,
                     CIcov.NDE.NC = cov.NDE.NC,
                     CIcov.NIE.NC = cov.NIE.NC,
                     CIcov.ATE.NC = cov.ATE.NC,
                     n = as.character(RECAP1[1,]$n),
                     pY1 = as.character(RECAP1[1,]$pY1),
                     A =  as.character(RECAP1[1,]$A),
                     beta_1 = beta_1,
                     ATE.true = RECAP1[1,]$ATE.true,
                     NDE.true = RECAP1[1,]$NDE.true,
                     NIE.true = RECAP1[1,]$NIE.true,
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
Res = as.data.frame(Res)
ColNames = colnames(Res[,c(1:38, 40:57)])
Res[ColNames] = sapply(Res[ColNames], as.numeric)
save(Res, file = paste0("Res_RandomizedUnblinded-beta1",
                        round(beta_1, digits = 3), ".RData"))
Res[ColNames] = round(Res[ColNames], digits = 3)
write.csv(Res, file = paste0("Res_RandomizedUnblinded-beta1",
                             round(beta_1, digits = 3), ".csv"))
write.xlsx(Res, paste0("Res_RandomizedUnblinded-beta1",
                       round(beta_1, digits = 3), ".xlsx"))

Res1 = Res[which(Res$n == 10000),]
latextable = cbind(Res1[, c(7:9,22,24,26,23,25,27,38,39,44)]) # only Joint-NC and the scenarios with n = 10,000
print(xtable(latextable, digits = 3), include.rownames=FALSE) # print the latex table
# Res2 = Res[which(Res$n == 5000),]
# latextable2 = cbind(Res2[, c(7:9,22,24,26,23,25,27,38,39,44)])
# print(xtable(latextable2, digits = 3), include.rownames=FALSE) # latex table for n = 5,000

## Plotting the results --------------------------------------------------------
RECAP$pY1 = factor(RECAP$pY1,
                   labels = c(expression(P(Y[1] == 1) == 0.025),
                              expression(P(Y[1] == 1) == 0.05),
                              expression(P(Y[1] == 1) == 0.14)))
RECAP$A = factor(RECAP$A,
                 labels = c(expression(a[low]==~0~","~a[medium]==~0.75~","~a[high]==~1.5),
                            expression(a[low]==~0~","~a[medium]==~1~","~a[high]==~2),
                            expression(a[low]==~0~","~a[medium]==~1~","~a[high]==~2.5) ))
RECAP1 = RECAP[which(RECAP$Approach == "NC.NIE" | RECAP$Approach == "NC.NDE" | RECAP$Approach == "NC.ATE"),]
RECAP1$Approach = factor(RECAP1$Approach)
RECAP1$Approach = factor(RECAP1$Approach,
                         labels = c(expression(hat(beta)^"*", hat(beta),hat(beta^"*"))))
plot = ggplot(RECAP1, aes(x = n, y = Effect.hat, color = Approach)) +
  geom_boxplot(coef = NULL) +
  geom_hline(aes(yintercept = NDE.true)) +
  geom_hline(aes(yintercept = ATE.true), linetype = "dotdash") +
  geom_hline(aes(yintercept = NIE.true), linetype = "dotted") + theme_light() +
  ylab("Estimate") + xlab("Unblinded randomized controlled trial size") +
  facet_grid(pY1~A, labeller = label_parsed, scales = "free_y") +
  theme(plot.title = element_text(size = 11),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 11),
        strip.background = element_rect(color="black", fill="white", size = 0.5, linetype="solid"),
        strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black")) +
  scale_colour_manual(values = c("deeppink", "cadetblue", "darkorange"),
                      labels = expression(hat(beta[1])^"*", hat(beta[1])^"*"-hat(beta[2])^"*",hat(beta[2])^"*"),
                      name = "Joint-NC estimate")
pdf(paste0("Comparison_RandomizedUnblinded-methodNC-beta1", round(beta_1, digits = 3), ".pdf"),  width = 10, height = 7)
plot
dev.off()
