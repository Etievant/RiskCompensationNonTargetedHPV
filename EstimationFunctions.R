## -----------------------------------------------------------------------------
## Description: This script contrains the functions to estimate the treatment 
##              effect and its variance, with the different methods proposed by 
##              Etievant et al. (Biometrics, 2023) and needed to replicate the 
##              simulation studies in Etievant (2025) 
##
##              Functions:
##
##  JointNC     Estimation with method Joint-NC. Relies on the joint estimation 
##              of the treatment effect on the targeted (i.e., primary) and 
##              untargeted (i.e., secondary) outcomes. Information on potential 
##              observed covariates or confounders is not used 
##
##  JointMH     Estimation with method Joint-MH. Relies on the joint estimation 
##              of the treatment effect on the targeted (i.e., primary) and 
##              untargeted (i.e., secondary) outcomes, using stratification on 
##              observed confounders with Mantel-Haenszel weights
##
##  JointReg    Estimation with method Joint-Reg. Relies on joint estimation of 
##              the treatment effect on the targeted (i.e., primary) and 
##              untargeted (i.e., secondary) outcomes, using regression models 
##              where observed confounders are included. For continuous 
##              confounders, quadratic polynomial functions are used
##
## JointRegCat  Estimation with method Joint-Reg, just as in function JointReg, 
##              but when observed confounders are all categorical. This function 
##              will be called by JointReg if argument Wcont is NULL
##
## JointRegCont Estimation with method Joint-Reg, just as in function JointReg, 
##              but when observed confounders are all continuous. This function 
##              will be called by JointReg if argument Wcat is NULL
## -----------------------------------------------------------------------------

### load package ---------------------------------------------------------------
library(gee)


### JointNC --------------------------------------------------------------------

## Arguments:
#   Y1  outcome of interest, assumed to be binary
#   Y2  secondary outcome, assumed to be categorical
#   T   treatment of interest, assumed to be binary. A log-link is used to 
#       relate Y1 and T and to relate Y2 and T

## Values:
#   beta_1.hat            final (``debiased'') estimate of the treatment effect 
#                         on Y1
#   se.beta_1.hat         sandwich standard error for beta_1.hat
#   beta_1.hat.naive      naive estimate of the treatment effect on Y1
#   beta_2.hat.naive      naive estimate of the treatment effect on Y2
#   se.beta_1.hat.naive   sandwich standard error for beta_1.hat.naive
#   se.beta_2.hat.naive   sandwich standard error for beta_2.hat.naive

JointNC = function(Y1 = Y1, Y2 = Y2, T = T){
  n                   = length(Y1) # sample size
  
  # Parameter estimation
  e.mu_1.tilde.hat    = sum((1 - T) * Y1) / sum(1 - T) # estimation of the intercept
  e.beta_1.tilde.hat  = sum(T * Y1) / sum(T) / e.mu_1.tilde.hat # naive estimation of the treatment effect on Y1
  e.mu_2.tilde.hat    = sum((1 - T) * Y2) / sum(1 - T) # estimation of the intercept
  e.beta_2.tilde.hat  = sum(T * Y2) / sum(T) / e.mu_2.tilde.hat # naive estimation of the treatment effect on Y2
  beta_1.tilde.hat    = log(e.beta_1.tilde.hat) 
  beta_2.tilde.hat    = log(e.beta_2.tilde.hat)
  # using the estimated non-zero treatment effect on the secondary outcome to reduce confounding bias in the estimated treatment effect on the primary outcome
  beta_1.hat          = beta_1.tilde.hat - beta_2.tilde.hat # ``de-biased'' treatment effect on Y1
  
  # Sandwich variance estimation 
  p1.tilde.hat      = e.mu_1.tilde.hat * (e.beta_1.tilde.hat * T + 1 - T)
  p2.tilde.hat      = e.mu_2.tilde.hat * (e.beta_2.tilde.hat * T + 1 - T)
  U.hat             = rbind((Y1 - p1.tilde.hat) / (1 - p1.tilde.hat), T * (Y1 - p1.tilde.hat) / (1 - p1.tilde.hat), Y2 - p2.tilde.hat, T * (Y2 - p2.tilde.hat))
  drond_U.hat       = matrix(NA, nrow = 4, ncol = 4)
  drond_U.hat[1,]   = rbind(sum(p1.tilde.hat * (Y1 - 1) / (1 - p1.tilde.hat)^2), sum(T * p1.tilde.hat * (Y1 - 1) / (1 - p1.tilde.hat)^2), 0, 0) / n
  drond_U.hat[2,]   = rbind(sum(T * p1.tilde.hat * (Y1 - 1) / (1 - p1.tilde.hat)^2), sum(T^2 * p1.tilde.hat * (Y1 - 1) / (1 - p1.tilde.hat)^2), 0, 0) / n
  drond_U.hat[3,]   = rbind(0, 0, - sum(p2.tilde.hat), - sum(T * p2.tilde.hat)) / n
  drond_U.hat[4,]   = rbind(0, 0, - sum(T * p2.tilde.hat), - sum(T * p2.tilde.hat)) / n
  if(is.na(e.beta_1.tilde.hat)){
    print("No untreated cases")
    beta_1.hat      = NA
    Var_theta.hat   = matrix(NA, nrow = 4, ncol = 4)
  }else{
    if(e.beta_1.tilde.hat == 0){
      print("No treated cases")
      beta_1.hat    = NA
      Var_theta.hat = matrix(NA, nrow = 4, ncol = 4)
    }else{
      if(e.beta_1.tilde.hat == Inf){
        print("No untreated cases")
        beta_1.hat    = NA
        Var_theta.hat = matrix(NA, nrow = 4, ncol = 4)
      }else{
        Var_theta.hat = 1 / n * solve(drond_U.hat) %*% (U.hat %*% t(U.hat) / n) %*% t(solve(drond_U.hat))
      }}}
  
  return(list(beta_1.hat = beta_1.hat, se.beta_1.hat = sqrt(Var_theta.hat[2,2] + Var_theta.hat[4,4] - 2 * Var_theta.hat[4,2]), beta_1.hat.naive = beta_1.tilde.hat, beta_2.hat.naive = beta_2.tilde.hat, se.beta_1.hat.naive = sqrt(Var_theta.hat[2,2]), se.beta_2.hat.naive = sqrt(Var_theta.hat[4,4])))
}
### ----------------------------------------------------------------------------


### JointMH --------------------------------------------------------------------

## Arguments:
#   Y1  outcome of interest, assumed to be binary
#   Y2  secondary outcome, assumed to be categorical
#   T   treatment of interest, assumed to be binary. A log-link is used to 
#       relate Y1 and T and to relate Y2 and T
#   W   categorical confounders used for adjustment

## Values:
#   beta_1.hat            final (``debiased'') estimate of the treatment effect 
#                         on Y1
#   se.beta_1.hat         sandwich standard error for beta_1.hat
#   beta_1.hat.naive      naive estimate of the treatment effect on Y1
#   beta_2.hat.naive      naive estimate of the treatment effect on Y2
#   se.beta_1.hat.naive   sandwich standard error for beta_1.hat.naive
#   se.beta_2.hat.naive   sandwich standard error for beta_2.hat.naive


JointMH = function(Y1 = Y1, Y2 = Y2, T = T, W = NULL){
  if(is.null(W)){
    # in the absence of measured covariates, perform the joint approach with no covariates (JointNC)
    JointNC(Y1 = Y1, Y2 = Y2, T = T)
  }else{
    n                 = length(Y1)
    W                 = as.factor(W)
    levelsW           = levels(W) # the strata
    K                 = length(levelsW)
    
    # Parameter estimation
    res               = NULL
    for(l in levelsW){
      nl              = sum(W == l) # number of individuals in the stratum
      nl1             = sum((W == l)&(T == 1)) # number of individuals who received the treatment in the stratum
      nl0             = sum((W == l)&(T == 0)) # number of individuals who did not receive the treatment in the stratum
      pl              = nl0 * nl1 / nl
      pl1             = nl0 / nl 
      pl0             = nl1 / nl
      Num_Y1          = sum(Y1[(W == l)&(T == 1)]) # number of treated cases in the stratum, for the primary outcome
      Denom_Y1        = sum(Y1[(W == l)&(T == 0)]) # number of untreated cases in the stratum, for the primary outcome
      Num_Y2          = sum(Y2[(W == l)&(T == 1)]) # number of treated cases in the stratum, for the secondary outcome
      Denom_Y2        = sum(Y2[(W == l)&(T == 0)]) # number of untreated cases in the stratum, for the secondary outcome
      res             = rbind(res, cbind(Num_Y1 = Num_Y1, Denom_Y1 = Denom_Y1, Num_Y2 = Num_Y2, Denom_Y2 = Denom_Y2, pl = pl, pl1 = pl1, pl0 = pl0))
    }
    res               = as.data.frame(res)
    if(sum(is.na(res$pl))){
      res             = res[-which(is.na(res$pl)), ] 
    }
    if(sum(res$pl == 0)){
      res             = res[-which(res$pl == 0), ] 
    }
    beta_1.tilde.hat  = log(sum(res$Num_Y1 * res$pl1) / sum(res$Denom_Y1 * res$pl0)) # estimated treatment effect on the primary outcome via stratification on W with MH weights
    beta_2.tilde.hat  = log(sum(res$Num_Y2 * res$pl1) / sum(res$Denom_Y2 * res$pl0)) # estimated treatment effect on the secondary outcome via stratification on W with MH weights
    # using the estimated non-zero treatment effect on the secondary outcome to reduce confounding bias in the estimated treatment effect on the primary outcome
    beta_1.hat        = beta_1.tilde.hat - beta_2.tilde.hat # ``de-biased'' treatment effect on Y1
    
    if(beta_1.tilde.hat == -Inf){
        print("No treated cases")
        beta_1.hat    = NA
        Var_theta.hat = matrix(NA, nrow = 2, ncol = 2)
    }else{
      if(beta_1.tilde.hat == Inf){
          print("No untreated cases")
          beta_1.hat    = NA
          Var_theta.hat = matrix(NA, nrow = 2, ncol = 2)
      }else{
          # Sandwich variance estimation
          U.hat             = rbind(res$pl1 * res$Num_Y1 - exp(beta_1.tilde.hat) * res$pl0 * res$Denom_Y1, res$pl1 * res$Num_Y2 - exp(beta_2.tilde.hat) * res$pl0 * res$Denom_Y2)
          bar_U.hat         = rowMeans(U.hat)
          drond_U.hat       = matrix(NA, nrow = 2, ncol = 2)
          drond_U.hat[1,]   = rbind(- exp(beta_1.tilde.hat) * sum(res$Denom_Y1 * res$pl0), 0) / K
          drond_U.hat[2,]   = rbind(0, - exp(beta_2.tilde.hat) * sum(res$Denom_Y2 * res$pl0)) / K
          Var_theta.hat     = 1 / (K) * solve(drond_U.hat) %*% ((U.hat - bar_U.hat) %*% t(U.hat - bar_U.hat) / (K - 1)) %*% t(solve(drond_U.hat))
        }}
        
    return(list(beta_1.hat = beta_1.hat, se.beta_1.hat = sqrt(Var_theta.hat[1,1] + Var_theta.hat[2,2] - 2*Var_theta.hat[1,2]), beta_1.hat.naive = beta_1.tilde.hat, beta_2.hat.naive = beta_2.tilde.hat, se.beta_1.hat.naive = sqrt(Var_theta.hat[1,1]), se.beta_2.hat.naive = sqrt(Var_theta.hat[2,2])))
  }
}
### ----------------------------------------------------------------------------


### JointReg -------------------------------------------------------------------

## Arguments:
#   Y1  outcome of interest, assumed to be binary
#   Y2  secondary outcome, assumed to be categorical
#   T   treatment of interest, assumed to be binary. A log-link is used to 
#       relate Y1 and T and to relate Y2 and T
# Wcat  categorical confounder used for adjustment
# Wcont continuous confounders used for adjustment. We assume Wcont affects Y1 
#       and Y2 through linear quadratic functions

## Values:
#   beta_1.hat            final (``debiased'') estimate of the treatment effect 
#                         on Y1
#   se.beta_1.hat         sandwich standard error for beta_1.hat
#   beta_1.hat.naive      naive estimate of the treatment effect on Y1
#   beta_2.hat.naive      naive estimate of the treatment effect on Y2
#   se.beta_1.hat.naive   sandwich standard error for beta_1.hat.naive
#   se.beta_2.hat.naive   sandwich standard error for beta_2.hat.naive


JointReg = function(Y1 = Y1, Y2 = Y2, T = T, Wcont = NULL, Wcat = NULL){
  n = length(Y1)
  
  if((is.null(Wcat))|(is.null(Wcont))){
    if((is.null(Wcat))&(is.null(Wcont))){
      # in the absence of measured covariates, perform the joint approach with no covariates (JointNC)
      JointNC(Y1 = Y1, Y2 = Y2, T = T)
    }else{
      if((is.null(Wcat))){
        JointRegCont(Y1 = Y1, Y2 = Y2, T = T, Wcont = Wcont)
      }else{
        JointRegCat(Y1 = Y1, Y2 = Y2, T = T, Wcat = Wcat)
      }}}
  else{
    Wcat                = as.matrix(Wcat)
    Wcont               = as.matrix(Wcont)
    Wcont2              = Wcont^2 # using a quadradic term for the continuous covariates
    levelsW             = levels(as.factor(Wcat)) # building dummy variables for the categories of the categorical covariates
    m                   = length(levelsW)
    WW                  = matrix(0, nrow = n, ncol = m)
    for(i in 1:m){
      WW[which(Wcat == levelsW[i]),i] = 1
    }
    if(sum((colSums(WW) == 0))){
      levelsWW          = levelsW[-which(colSums(WW) == 0)]
      WW                = WW[,-which(colSums(WW) == 0)]
      m                 = ncol(WW)
    }
    WW                  = WW[,-m]
    WW                  = as.matrix(WW)
    
    # regression models for the primary and secondary outcome
    # alternatively, one could have solved the associated estimating equations
    mod1.gee            = glm(Y1 ~ T + Wcont + Wcont2 + WW, family = binomial(link = "log")) # binomial log link as Y1 is binary
    mod2.gee            = glm(Y2 ~ T + Wcont + Wcont2 + WW, family = poisson) # poisson link as Y2 is categorical
    
    # Parameter estimation
    mu_1.hat            = mod1.gee$coefficients[1]
    beta_1.hat          = mod1.gee$coefficients[2]
    alpha_1cont.hat     = c(mod1.gee$coefficients[3:(2 + ncol(Wcont) * 2)])
    alpha_1cat.hat      = c(mod1.gee$coefficients[(3 + ncol(Wcont) * 2):(length(mod1.gee$coefficients))])
    mu_2.hat            = mod2.gee$coefficients[1]
    beta_2.hat          = mod2.gee$coefficients[2]
    alpha_2cont.hat     = c(mod2.gee$coefficients[3:(2 + ncol(Wcont) * 2)])
    alpha_2cat.hat      = c(mod2.gee$coefficients[(3 + ncol(Wcont) * 2):(length(mod2.gee$coefficients))])
    
    # Sandwich variance estimation
    p1.hat              = exp(mu_1.hat + beta_1.hat * T + cbind(Wcont,Wcont2) %*% alpha_1cont.hat + WW %*% alpha_1cat.hat) 
    p2.hat              = exp(mu_2.hat + beta_2.hat * T + cbind(Wcont,Wcont2) %*% alpha_2cont.hat + WW %*% alpha_2cat.hat) 
    
    ee1                       = (Y1 - p1.hat) / (1 - p1.hat)             # first estimating equation for Y1 (i.e., for the intercept)
    ee1T                      = T * (Y1 - p1.hat) / (1 - p1.hat)         # "T" estimating equation for Y1
    ee1Wcont                  = cbind(Wcont, Wcont2) * matrix(rep((Y1 - p1.hat) / (1 - p1.hat), ncol(Wcont) * 2 ), nrow = n, byrow = F) # "Wcont" estimating equation for Y1
    ee1wcat                   = as.matrix((Y1 - p1.hat) / (1 - p1.hat))  # used for the "Wcat" estimating equation for Y1
    deriv.ee1TWcont           = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # derivative of the "T" estimating equation (for Y1) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to T
    deriv.ee1Wcontwcat        = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # used for the derivative of the "Wcont" estimating equation (for Y1) wr to Wcat, or the derivative of the "Wcat" estimating equations wr to "Wcont"
    deriv.ee1Twcat            = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2) # used for the derivative of the "T" estimating equation (for Y1) wr to Wcat, or the derivative of the "Wcat" estimating equations wr to T. Will be also used for the derivative of the "Wcat" estimating equation wr to "Wcat"
    deriv.ee1wcat             = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # used for the derivative of the first estimating equation (for Y1) wr to "Wcat", or the derivative of the "Wcat" estimating equations wr to the intercept
    deriv.ee1wcont            = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # derivative of the first estimating equation (for Y1) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to the intercept. Will be also used for the derivative of the "Wcont" estimating equation wr to "Wcont"
    deriv.ee1T                = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2) # derivative of the first estimating equation (for Y1) wr to T, or the derivative of the "T" estimating equations wr to the intercept
    deriv.ee1                 = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # derivative of the first estimating equation (for Y1) wr to the intercept
    
    ee2                       = (Y2 - p2.hat) # same for Y2
    ee2T                      = T * (Y2 - p2.hat)
    ee2Wcont                  = cbind(Wcont, Wcont2) * matrix(rep((Y2 - p2.hat), ncol(Wcont) * 2 ), nrow = n, byrow = F) 
    ee2wcat                   = as.matrix((Y2 - p2.hat))            
    deriv.ee2TWcont           = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * T * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F) 
    deriv.ee2Wcontwcat        = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F) 
    deriv.ee2Twcat            = as.matrix(p2.hat * T * (- 1))         
    deriv.ee2wcat             = as.matrix(p2.hat * (- 1))            
    deriv.ee2wcont            = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F)
    deriv.ee2T                = as.matrix(p2.hat * T * (- 1))
    deriv.ee2                 = as.matrix(p2.hat * (- 1))
    
    ee1Wcat                   = WW    # "Wcat" estimating equation for Y1
    ee2Wcat                   = WW    # same for Y2
    deriv.ee1Wcat             = WW    # derivative of the T estimating equation (for Y1) wr to Wcat, or derivative of Wcat estimating equations wr to T coefficient
    deriv.ee2Wcat             = WW    # same for Y2
    deriv.ee1TWcat            = WW    # derivative of the first estimating equation (for Y1) wr to Wcat, or derivative of the Wcat estimating equations wr to the intercept
    deriv.ee2TWcat            = WW    # same for Y2
    deriv.ee1WcontWcat        = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(WW))    # derivative of the Wcat estimating equation (for Y1) wr to Wcont, or derivative of the Wcont estimating equations wr to Wcat
    deriv.ee2WcontWcat        = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(WW))    # same for Y2
    for(i in 1:(m-1)){
      ee1Wcat[,i]             = WW[,i] * ee1wcat 
      ee2Wcat[,i]             = WW[,i] * ee2wcat 
      deriv.ee1Wcat[,i]       = WW[,i] * deriv.ee1wcat 
      deriv.ee2Wcat[,i]       = WW[,i] * deriv.ee2wcat 
      deriv.ee1TWcat[,i]      = WW[,i] * deriv.ee1Twcat 
      deriv.ee2TWcat[,i]      = WW[,i] * deriv.ee2Twcat 
      deriv.ee1WcontWcat[,i]  = colSums(WW[,i] * deriv.ee1Wcontwcat)
      deriv.ee2WcontWcat[,i]  = colSums(WW[,i] * deriv.ee2Wcontwcat)
    }
    deriv.ee1Wcont        = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(Wcont) * 2) # derivative of the "Wcont" estimating equation for Y1
    deriv.ee2Wcont        = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(Wcont) * 2) # same for Y2
    for(i in 1:(ncol(Wcont) * 2)){
      deriv.ee1Wcont[,i]  = colSums(cbind(Wcont,Wcont2)[,i] * deriv.ee1wcont)
      deriv.ee2Wcont[,i]  = colSums(cbind(Wcont,Wcont2)[,i] * deriv.ee2wcont)
    }
    U.hat                 = cbind(ee1, ee1T, ee1Wcont, ee1Wcat, ee2, ee2T, ee2Wcont, ee2Wcat) # estimating function evaluated at theta.hat for each individual
    U.hat.mean            = matrix(rep(colMeans(U.hat), n), nrow = n, byrow = TRUE)
    
    # Derivative of the estimating equation evaluated in theta.hat
    # used for the sandwich variance
    drond_U.hat                           = matrix(NA, nrow = ncol(U.hat), ncol = ncol(U.hat))
    drond_U.hat[1,]                       = c(sum(deriv.ee1), sum(deriv.ee1T), colSums(deriv.ee1wcont), colSums(deriv.ee1Wcat), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[2,]                       = c(sum(deriv.ee1T), sum(deriv.ee1T), colSums(deriv.ee1TWcont), colSums(deriv.ee1TWcat), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[3:(2 + ncol(Wcont) * 2),] = cbind(colSums(deriv.ee1wcont), colSums(deriv.ee1TWcont), deriv.ee1Wcont, deriv.ee1WcontWcat, matrix(rep(0, ncol(Wcont) * 2 * ncol(U.hat) / 2), nrow = ncol(Wcont) * 2, byrow = F)) / n 
    
    UU1                                                       = matrix(0, nrow = (m-1), ncol = 2*(2 +  ncol(Wcont) * 2 + m-1))
    UU1[,1]                                                   = colSums(deriv.ee1Wcat) / n
    UU1[,2]                                                   = colSums(deriv.ee1TWcat) / n
    UU1[,3:(2 + ncol(Wcont) * 2)]                             = t(deriv.ee1WcontWcat) / n
    UU1[,(3 + ncol(Wcont) * 2):(2 +  ncol(Wcont) * 2 + m-1)]  = diag(colSums(deriv.ee1Wcat), nrow = length((3 + ncol(Wcont) * 2):(2 +  ncol(Wcont) * 2 + m-1))) / n
    drond_U.hat[(3 + ncol(Wcont) * 2):(2 +  ncol(Wcont) * 2 + m-1),] = UU1
    
    drond_U.hat[1 + ncol(U.hat) / 2,]                                           = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2), sum(deriv.ee2T), colSums(deriv.ee2wcont), colSums(deriv.ee2Wcat)) / n
    drond_U.hat[2 + ncol(U.hat) / 2,]                                           = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2T), sum(deriv.ee2T), colSums(deriv.ee2TWcont), colSums(deriv.ee2TWcat)) / n
    drond_U.hat[(3 + ncol(U.hat) / 2):(2 + ncol(Wcont) * 2 + ncol(U.hat) / 2),] = cbind(matrix(rep(0, ncol(Wcont) * 2 * ncol(U.hat) / 2), nrow = ncol(Wcont) * 2, byrow = F), colSums(deriv.ee2wcont), colSums(deriv.ee2TWcont), deriv.ee2Wcont, deriv.ee2WcontWcat) / n 
    
    UU2                                                                                                   = matrix(0, nrow = (m-1), ncol = 2*(2 +  ncol(Wcont) * 2 + m-1))
    UU2[,ncol(U.hat) / 2 + 1]                                                                             = colSums(deriv.ee2Wcat) / n
    UU2[,ncol(U.hat) / 2 + 2]                                                                             = colSums(deriv.ee2TWcat) / n
    UU2[,(3 + ncol(U.hat) / 2):(2 + ncol(U.hat) / 2 + ncol(Wcont) * 2)]                                   = t(deriv.ee2WcontWcat) / n
    UU2[,(3 + ncol(U.hat) / 2 + ncol(Wcont) * 2):(2 + ncol(U.hat) / 2 +  ncol(Wcont) * 2 + m-1)]          = diag(colSums(deriv.ee2Wcat), nrow = length((3 + ncol(U.hat) / 2 + ncol(Wcont) * 2):(2+ ncol(U.hat) / 2 +  ncol(Wcont) * 2 + m-1))) / n
    drond_U.hat[(3 + ncol(U.hat) / 2 + ncol(Wcont) * 2):(2 + ncol(U.hat) / 2 +  ncol(Wcont) * 2 + m-1),]  = UU2
    
    # Sandwich variance 
    Var_theta.hat = 1 / n * solve(drond_U.hat) %*% (t(U.hat - U.hat.mean) %*% (U.hat - U.hat.mean) / n) %*% t(solve(drond_U.hat)) 
  
    return(list(beta_1.hat = beta_1.hat - beta_2.hat, se.beta_1.hat = sqrt(Var_theta.hat[2,2] + Var_theta.hat[(2 + ncol(U.hat) / 2),(2 + ncol(U.hat) / 2)] - 2 * Var_theta.hat[(2 + ncol(U.hat) / 2),2]), beta_1.hat.naive = beta_1.hat, beta_2.hat.naive = beta_2.hat, se.beta_1.hat.naive = sqrt(Var_theta.hat[2,2]), se.beta_2.hat.naive = sqrt(Var_theta.hat[(2 + ncol(U.hat) / 2),(2 + ncol(U.hat) / 2)])))
  }
}
### ----------------------------------------------------------------------------


### JointRegCont ---------------------------------------------------------------

## Arguments:
#   Y1  outcome of interest, assumed to be binary
#   Y2  secondary outcome, assumed to be categorical
#   T   treatment of interest, assumed to be binary. A log-link is used to 
#       relate Y1 and T and to relate Y2 and T
# Wcont continuous confounders used for adjustment. We assume Wcont affects Y1 
#       and Y2 through linear quadratic functions

## Values:
#   beta_1.hat            final (``debiased'') estimate of the treatment effect 
#                         on Y1
#   se.beta_1.hat         sandwich standard error for beta_1.hat
#   beta_1.hat.naive      naive estimate of the treatment effect on Y1
#   beta_2.hat.naive      naive estimate of the treatment effect on Y2
#   se.beta_1.hat.naive   sandwich standard error for beta_1.hat.naive
#   se.beta_2.hat.naive   sandwich standard error for beta_2.hat.naive


JointRegCont = function(Y1 = Y1, Y2 = Y2, T = T, Wcont = NULL){
  n                   = length(Y1)
  
  if(is.null(Wcont)){
    JointNC(Y1 = Y1, Y2 = Y2, T = T)
  }else{
    
    Wcont             = as.matrix(Wcont)
    Wcont2            = Wcont^2 # using a quadradic term for the continuous covariates
    
    # regression models for the primary and secondary outcome
    # alternatively, one could have solved the associated estimating equations
    mod1.gee          = glm(Y1 ~ T + Wcont + Wcont2 , family = binomial(link = "log"))
    mod2.gee          = glm(Y2 ~ T + Wcont + Wcont2, family = poisson)
    
    # Parameter estimation
    mu_1.hat          = mod1.gee$coefficients[1]
    beta_1.hat        = mod1.gee$coefficients[2]
    alpha_1cont.hat   = c(mod1.gee$coefficients[3:(2 + ncol(Wcont) * 2)])
    mu_2.hat          = mod2.gee$coefficients[1]
    beta_2.hat        = mod2.gee$coefficients[2]
    alpha_2cont.hat   = c(mod2.gee$coefficients[3:(2 + ncol(Wcont) * 2)])
    
    # Sandwich variance estimation
    p1.hat          = exp(mu_1.hat + beta_1.hat * T + cbind(Wcont,Wcont2) %*% alpha_1cont.hat) 
    p2.hat          = exp(mu_2.hat + beta_2.hat * T + cbind(Wcont,Wcont2) %*% alpha_2cont.hat) 
    
    ee1                     = (Y1 - p1.hat) / (1 - p1.hat)             # first estimating equation for Y1
    ee1T                    = T * (Y1 - p1.hat) / (1 - p1.hat)         # "T" estimating equation for Y1
    ee1Wcont                = cbind(Wcont, Wcont2) * matrix(rep((Y1 - p1.hat) / (1 - p1.hat), ncol(Wcont) * 2 ), nrow = n, byrow = F)         # "Wcont" estimating equation for Y1
    deriv.ee1TWcont         = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # derivative of the "T" estimating equation (for Y1) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to T
    deriv.ee1wcont          = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # derivative of the first estimating equation (for Y1) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to the intercept. Will be also used for the derivative of the "Wcont" estimating equation wr to "Wcont"
    deriv.ee1T              = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2)  # derivative of the first estimating equation (for Y1) wr to T, or the derivative of the "T" estimating equations wr to the intercept
    deriv.ee1               = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # derivative of the first estimating equation (for Y1) wr to the intercept.
    
    ee2                     = (Y2 - p2.hat)                            # same for Y2
    ee2T                    = T * (Y2 - p2.hat)                   
    ee2Wcont                = cbind(Wcont, Wcont2) * matrix(rep((Y2 - p2.hat), ncol(Wcont) * 2 ), nrow = n, byrow = F) 
    deriv.ee2TWcont         = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * T * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F)
    deriv.ee2wcont          = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F) 
    deriv.ee2T              = as.matrix(p2.hat * T * (- 1)) 
    deriv.ee2               = as.matrix(p2.hat * (- 1))
    
    deriv.ee1Wcont = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(Wcont) * 2) # derivative of the "Wcont" estimating equation for Y1
    deriv.ee2Wcont = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(Wcont) * 2) # same for Y2
    for(i in 1:(ncol(Wcont) * 2)){
      deriv.ee1Wcont[,i] = colSums(cbind(Wcont,Wcont2)[,i] * deriv.ee1wcont)
      deriv.ee2Wcont[,i] = colSums(cbind(Wcont,Wcont2)[,i] * deriv.ee2wcont)
    }
    
    U.hat             = cbind(ee1, ee1T, ee1Wcont, ee2, ee2T, ee2Wcont) # estimating function evaluated in theta.hat for each individual
    U.hat.mean        = matrix(rep(colMeans(U.hat), n), nrow = n, byrow = TRUE)
    
    # derivative of the estimating equation evaluated in theta.hat
    # used for the sandwich variance
    drond_U.hat       = matrix(NA, nrow = ncol(U.hat), ncol = ncol(U.hat))
    drond_U.hat[1,]   = c(sum(deriv.ee1), sum(deriv.ee1T), colSums(deriv.ee1wcont), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[2,]   = c(sum(deriv.ee1T), sum(deriv.ee1T), colSums(deriv.ee1TWcont), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[3:(2 + ncol(Wcont) * 2),]                                           = cbind(colSums(deriv.ee1wcont), colSums(deriv.ee1TWcont), deriv.ee1Wcont,  matrix(rep(0, ncol(Wcont) * 2 * ncol(U.hat) / 2), nrow = ncol(Wcont) * 2, byrow = F)) / n 
    drond_U.hat[1 + ncol(U.hat) / 2,]                                               = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2), sum(deriv.ee2T), colSums(deriv.ee2wcont)) / n
    drond_U.hat[2 + ncol(U.hat) / 2,]                                               = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2T), sum(deriv.ee2T), colSums(deriv.ee2TWcont)) / n
    drond_U.hat[(3 + ncol(U.hat) / 2):(2 + ncol(Wcont) * 2 + ncol(U.hat) / 2),]     = cbind(matrix(rep(0, ncol(Wcont) * 2 * ncol(U.hat) / 2), nrow = ncol(Wcont) * 2, byrow = F), colSums(deriv.ee2wcont), colSums(deriv.ee2TWcont), deriv.ee2Wcont) / n 
    
    # Sandwich variance
    Var_theta.hat = 1 / n * solve(drond_U.hat) %*% (t(U.hat - U.hat.mean) %*% (U.hat - U.hat.mean) / n) %*% t(solve(drond_U.hat)) 
    
    return(list(beta_1.hat = beta_1.hat - beta_2.hat, se.beta_1.hat = sqrt(Var_theta.hat[2,2] + Var_theta.hat[(2 + ncol(U.hat) / 2),(2 + ncol(U.hat) / 2)] - 2 * Var_theta.hat[(2 + ncol(U.hat) / 2),2]), beta_1.hat.naive = beta_1.hat, beta_2.hat.naive = beta_2.hat, se.beta_1.hat.naive = sqrt(Var_theta.hat[2,2]), se.beta_2.hat.naive = sqrt(Var_theta.hat[(2 + ncol(U.hat) / 2),(2 + ncol(U.hat) / 2)])))
  }
}
### ----------------------------------------------------------------------------


### JointRegCat ----------------------------------------------------------------

## Arguments:
#   Y1  outcome of interest, assumed to be binary
#   Y2  secondary outcome, assumed to be categorical
#   T   treatment of interest, assumed to be binary. A log-link is used to 
#       relate Y1 and T and to relate Y2 and T
# Wcat  categorical confounder used for adjustment

## Values:
#   beta_1.hat            final (``debiased'') estimate of the treatment effect 
#                         on Y1
#   se.beta_1.hat         sandwich standard error for beta_1.hat
#   beta_1.hat.naive      naive estimate of the treatment effect on Y1
#   beta_2.hat.naive      naive estimate of the treatment effect on Y2
#   se.beta_1.hat.naive   sandwich standard error for beta_1.hat.naive
#   se.beta_2.hat.naive   sandwich standard error for beta_2.hat.naive


JointRegCat = function(Y1 = Y1, Y2 = Y2, T = T, Wcat = NULL){
  
  n = length(Y1)
  
  if(is.null(Wcat)){
    JointNC(Y1 = Y1, Y2 = Y2, T = T)
  }else{
    
    Wcat          = as.matrix(Wcat)
    
    levelsW       = levels(as.factor(Wcat)) # building dummy variables for the categories of the categorical covariates
    m             = length(levelsW)
    WW            = matrix(0, nrow = n, ncol = m)
    for(i in 1:m){
      WW[which(Wcat == levelsW[i]),i] = 1
    }
    if(sum((colSums(WW) == 0))){
      levelsWW    = levelsW[-which(colSums(WW) == 0)]
      WW          = WW[,-which(colSums(WW) == 0)]
      m           = ncol(WW)
    }
    WW            = WW[,-m]
    WW            = as.matrix(WW)
    
    # regression models for the primary and secondary outcome
    # alternatively, one could have solved the associated estimating equations
    mod1.gee      = glm(Y1 ~ T + WW, family = binomial(link = "log"))
    mod2.gee      = glm(Y2 ~ T + WW, family = poisson)
    
    # Parameter estimation
    mu_1.hat      = mod1.gee$coefficients[1]
    beta_1.hat    = mod1.gee$coefficients[2]
    alpha_1cat.hat= c(mod1.gee$coefficients[(3):(length(mod1.gee$coefficients))])
    mu_2.hat      = mod2.gee$coefficients[1]
    beta_2.hat    = mod2.gee$coefficients[2]
    alpha_2cat.hat= c(mod2.gee$coefficients[(3):(length(mod2.gee$coefficients))])
    
    # Sandwich variance estimation
    p1.hat        = exp(mu_1.hat + beta_1.hat * T + WW %*% alpha_1cat.hat) 
    p2.hat        = exp(mu_2.hat + beta_2.hat * T + WW %*% alpha_2cat.hat) 
    
    ee1                     = (Y1 - p1.hat) / (1 - p1.hat)             # first estimating equation for Y1
    ee1T                    = T * (Y1 - p1.hat) / (1 - p1.hat)         # "T" estimating equation for Y1
    ee1wcat                 = as.matrix((Y1 - p1.hat) / (1 - p1.hat))  # used for the "Wcat" estimating equation for Y1
    deriv.ee1Twcat          = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2) # used for the derivative of the "T" estimating equation (for Y1) wr to Wcat, or the derivative of the "Wcat" estimating equations wr to T.  Will be also used for the derivative of the "Wcat" estimating equation wr to "Wcat"
    deriv.ee1wcat           = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # used for the derivative of the first estimating equation (for Y1) wr to "Wcat", or the derivative of the "Wcat" estimating equations wr to the intercept
    deriv.ee1T              = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2) # derivative of the first estimating equation (for Y1) wr to T, or the derivative of the "T" estimating equations wr to the intercept
    deriv.ee1               = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # derivative of the first estimating equation (for Y1) wr to the intercept.
    
    ee2                     = (Y2 - p2.hat)                            # same for Y2
    ee2T                    = T * (Y2 - p2.hat)
    ee2wcat                 = as.matrix((Y2 - p2.hat))  
    deriv.ee2Twcat          = as.matrix(p2.hat * T * (- 1))  
    deriv.ee2wcat           = as.matrix(p2.hat * (- 1)) 
    deriv.ee2T              = as.matrix(p2.hat * T * (- 1)) 
    deriv.ee2               = as.matrix(p2.hat * (- 1)) 
    
    ee1Wcat                 = WW    # "Wcat" estimating equation for Y1
    ee2Wcat                 = WW    # same for Y2
    deriv.ee1Wcat           = WW    # derivative of the T estimating equation (for Y1) wr to Wcat, or derivative of Wcat estimating equations wr to T coefficient.
    deriv.ee2Wcat           = WW    # same for Y2
    deriv.ee1TWcat          = WW    # derivative of the first estimating equation (for Y1) wr to Wcat, or derivative of the Wcat estimating equations wr to the intercept.
    deriv.ee2TWcat          = WW    # same for Y2
    for(i in 1:(m-1)){
      ee1Wcat[,i]           = WW[,i] * ee1wcat 
      ee2Wcat[,i]           = WW[,i] * ee2wcat 
      deriv.ee1Wcat[,i]     = WW[,i] * deriv.ee1wcat 
      deriv.ee2Wcat[,i]     = WW[,i] * deriv.ee2wcat 
      deriv.ee1TWcat[,i]    = WW[,i] * deriv.ee1Twcat 
      deriv.ee2TWcat[,i]    = WW[,i] * deriv.ee2Twcat 
    }
    U.hat             = cbind(ee1, ee1T, ee1Wcat, ee2, ee2T, ee2Wcat) # estimating function evaluated in theta.hat for each individual
    U.hat.mean        = matrix(rep(colMeans(U.hat), n), nrow = n, byrow = TRUE)
    
    # derivative of the estimating equation evaluated in theta.hat
    # used for the sandwich variance
    drond_U.hat       = matrix(NA, nrow = ncol(U.hat), ncol = ncol(U.hat))
    drond_U.hat[1,]   = c(sum(deriv.ee1), sum(deriv.ee1T), colSums(deriv.ee1Wcat), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[2,]   = c(sum(deriv.ee1T), sum(deriv.ee1T), colSums(deriv.ee1TWcat), rep(0, ncol(U.hat) / 2)) / n
    UU1               = matrix(0, nrow = (m-1), ncol = 2*(2 +  + m-1))
    UU1[,1]           = colSums(deriv.ee1Wcat) / n
    UU1[,2]           = colSums(deriv.ee1TWcat) / n
    UU1[,(3 ):(2  + m-1)]                                               = diag(colSums(deriv.ee1Wcat), nrow = length((3 ):(2  + m-1)), ncol = length((3 ):(2  + m-1))) / n
    drond_U.hat[(3 ):(2 + m-1),]                                        = UU1
    
    drond_U.hat[1 + ncol(U.hat) / 2,]                                   = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2), sum(deriv.ee2T), colSums(deriv.ee2Wcat)) / n
    drond_U.hat[2 + ncol(U.hat) / 2,]                                   = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2T), sum(deriv.ee2T),  colSums(deriv.ee2TWcat)) / n
    UU2                                                                 = matrix(0, nrow = (m-1), ncol = 2*(2 + m-1))
    UU2[,ncol(U.hat) / 2 + 1]                                           = colSums(deriv.ee2Wcat) / n
    UU2[,ncol(U.hat) / 2 + 2]                                           = colSums(deriv.ee2TWcat) / n
    UU2[,(3 + ncol(U.hat) / 2 ):(2+ ncol(U.hat) / 2  + m-1)]            = diag(colSums(deriv.ee2Wcat), nrow = length((3 + ncol(U.hat) / 2 ):(2+ ncol(U.hat) / 2  + m-1))) / n
    drond_U.hat[(3 + ncol(U.hat) / 2 ):(2 + ncol(U.hat) / 2 + m-1),]    = UU2
    
    # Sandwich variance
    Var_theta.hat = 1 / n * solve(drond_U.hat) %*% (t(U.hat - U.hat.mean) %*% (U.hat - U.hat.mean) / n) %*% t(solve(drond_U.hat)) 
    
    return(list(beta_1.hat = beta_1.hat - beta_2.hat, se.beta_1.hat = sqrt(Var_theta.hat[2,2] + Var_theta.hat[(2 + ncol(U.hat) / 2),(2 + ncol(U.hat) / 2)] - 2 * Var_theta.hat[(2 + ncol(U.hat) / 2),2]), beta_1.hat.naive = beta_1.hat, beta_2.hat.naive = beta_2.hat, se.beta_1.hat.naive = sqrt(Var_theta.hat[2,2]), se.beta_2.hat.naive = sqrt(Var_theta.hat[(2 + ncol(U.hat) / 2),(2 + ncol(U.hat) / 2)])))
  }
}
### ----------------------------------------------------------------------------
