## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=11.5, fig.height=8.5 
)

## ----sim_ctns0-----------------------------------------------------------
library(StratifiedMedicine)
dat_ctns = generate_subgrp_data(family="gaussian")
Y = dat_ctns$Y
X = dat_ctns$X # 50 covariates, 46 are noise variables, X1 and X2 are truly predictive
A = dat_ctns$A # binary treatment, 1:1 randomized 

## ----table_steps, echo=FALSE---------------------------------------------
library(knitr)
summ.table = data.frame( `Model` = c("filter", "ple", "submod", "param"),
                        `Required Outputs` = c("filter.vars", "list(mod,pred.fun)",
                                               "list(mod,pred.fun)", 
                                               "param.dat"),
                        `Description` = c("Variables that pass filter",
                                          "Model fit(s) and prediction function",
                                          "Model fit(s) and prediction function",
                                          "Parameter Estimates (overall and subgroups)") )
kable( summ.table, caption = "Key Outputs by Model" )

## ----user_filter_template------------------------------------------------
model_template(type="filter")

## ----user_filter---------------------------------------------------------
filter_lasso = function(Y, A, X, lambda="lambda.min", family="gaussian", ...){
  require(glmnet)
  ## Model matrix X matrix #
  X = model.matrix(~. -1, data = X )

  ##### Elastic Net on estimated ITEs #####
  set.seed(6134)
  if (family=="survival") { family = "cox"  }
  mod <- cv.glmnet(x = X, y = Y, nlambda = 100, alpha=1, family=family)

  ### Extract filtered variable based on lambda ###
  VI <- coef(mod, s = lambda)[,1]
  VI = VI[-1]
  filter.vars = names(VI[VI!=0])
  return( list(filter.vars=filter.vars) )
}

## ----user_ple_template---------------------------------------------------
model_template(type="ple")

## ----user_ple------------------------------------------------------------
ple_ranger_mtry = function(Y, A, X, Xtest, mtry=5, ...){
   require(ranger)
   ## Split data by treatment ###
    train0 =  data.frame(Y=Y[A==0], X[A==0,])
    train1 =  data.frame(Y=Y[A==1], X[A==1,])
    # Trt 0 #
    mod0 <- ranger(Y ~ ., data = train0, seed=1, mtry = mtry)
    # Trt 1 #
    mod1 <- ranger(Y ~ ., data = train1, seed=2, mtry = mtry)
    mod = list(mod0=mod0, mod1=mod1)
    pred.fun <- function(mod, X){
      mu1_hat <- predict( mod$mod1, X )$predictions
      mu0_hat <- predict( mod$mod0, X )$predictions
      mu_hat <- data.frame(mu1 = mu1_hat, mu0 = mu0_hat, PLE = mu1_hat-mu0_hat)
      return(mu_hat)
      }
    res = list(mod=mod, pred.fun=pred.fun)
    return( res )
}

## ----user_submod_template------------------------------------------------
model_template(type="submod")

## ----user_submod---------------------------------------------------------
submod_lmtree_pred = function(Y, A, X, Xtest, mu_train, ...){
  require(partykit)
  ## Fit Model ##
  mod <- lmtree(Y~A | ., data = X, parm=2) ##parm=2 focuses on treatment interaction #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
     Subgrps <- NULL
     Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
     return( list(Subgrps=Subgrps) )
  }
  ## Return Results ##
  return(  list(mod=mod, pred.fun=pred.fun) )
}

## ----user_param_template-------------------------------------------------
model_template(type="param")

## ----user_param----------------------------------------------------------

### Robust linear Regression: E(Y|A=1) - E(Y|A=0) ###
param_rlm = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, ...){
  require(MASS)
  indata = data.frame(Y=Y,A=A, X)
  mod.ovrl = rlm(Y ~ A , data=indata)
  param.dat0 = data.frame( Subgrps=0, N = dim(indata)[1],
                           estimand = "E(Y|A=1)-E(Y|A=0)",
                           est = summary(mod.ovrl)$coefficients[2,1],
                           SE = summary(mod.ovrl)$coefficients[2,2] )
  param.dat0$LCL = with(param.dat0, est-qt(1-alpha_ovrl/2, N-1)*SE)
  param.dat0$UCL = with(param.dat0, est+qt(1-alpha_ovrl/2, N-1)*SE)
  param.dat0$pval = with(param.dat0, 2*pt(-abs(est/SE), df=N-1) )

  ## Subgroup Specific Estimate ##
  looper = function(s){
    rlm.mod = tryCatch( rlm(Y ~ A , data=indata[Subgrps==s,]),
                       error = function(e) "param error" )
    n.s = dim(indata[Subgrps==s,])[1]
    est = summary(rlm.mod)$coefficients[2,1]
    SE = summary(rlm.mod)$coefficients[2,2]
    LCL =  est-qt(1-alpha_ovrl/2, n.s-1)*SE
    UCL =  est+qt(1-alpha_ovrl/2, n.s-1)*SE
    pval = 2*pt(-abs(est/SE), df=n.s-1)
    return( c(est, SE, LCL, UCL, pval) )
  }
  S_levels = as.numeric( names(table(Subgrps)) )
  S_N = as.numeric( table(Subgrps) )
  param.dat = lapply(S_levels, looper)
  param.dat = do.call(rbind, param.dat)
  param.dat = data.frame( S = S_levels, N=S_N, estimand="E(Y|A=1)-E(Y|A=0)", param.dat)
  colnames(param.dat) = c("Subgrps", "N", "estimand", "est", "SE", "LCL", "UCL", "pval")
  param.dat = rbind( param.dat0, param.dat)
  return( param.dat )
}


## ----user_prism, warnings=FALSE, message=FALSE---------------------------

res_user1 = PRISM(Y=Y, A=A, X=X, family="gaussian", filter="filter_lasso", 
             ple = "ple_ranger_mtry", submod = "submod_lmtree_pred",
             param="param_rlm")
## variables that remain after filtering ##
res_user1$filter.vars
## Subgroup model: lmtree searching for predictive only ##
plot(res_user1$submod.fit$mod)
## Parameter estimates/inference
res_user1$param.dat
plot(res_user1, type="forest")


