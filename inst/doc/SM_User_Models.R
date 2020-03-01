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
filter_template = function(Y, A, X, ...){
  # Step 1: Fit Filter Model #
  mod <- # model call 
  # Step 2: Extract variables that pass the filter #
  filter.vars <- # depends on mod fit
  # Return model fit and filtered variables #
  res = list(mod=mod, filter.vars=filter.vars)
  return( res )
}

## ----user_filter---------------------------------------------------------
filter_lasso = function(Y, A, X, lambda="lambda.min", family="gaussian", ...){
  require(glmnet)
  ## Model matrix X matrix #
  X = model.matrix(~. -1, data = X )

  ##### Elastic Net ##
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
ple_template <- function(Y, A, X, Xtest, ...){
  # Step 1: Fit PLE Model #
  # for example: Estimate E(Y|A=1,X), E(Y|A=0,X), E(Y|A=1,X)-E(Y|A=0,X)
  mod <- # ple model call 
  # mod = list(mod0=mod0, mod1=mod1) # If multiple fitted models, combine into list
  # Step 2: Predictions
  # Option 1: Create a Prediction Function #
  pred.fun <- function(mod, X){
    mu_hat <- # data-frame of predictions 
    return(mu_hat)
  }
  # Option 2: Directly Output Predictions (here, we still use pred.fun) #
  mu_train <- pred.fun(mod, X)
  mu_test <- pred.fun(mod, Xtest)
      
  # Return model fits and pred.fun (or just mu_train/mu_test) #
  res <- list(mod=mod, pred.fun=pred.fun, mu_train=mu_train, mu_test=mu_test)
  return( res )
}

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
      mu_1 <- predict( mod$mod1, X )$predictions
      mu_0 <- predict( mod$mod0, X )$predictions
      mu_hat <- data.frame(mu_1 = mu_1, mu_0 = mu_0, PLE = mu_1-mu_0)
      return(mu_hat)
      }
    res = list(mod=mod, pred.fun=pred.fun)
    return( res )
}

## ----user_submod_template------------------------------------------------
submod_template <- function(Y, A, X, Xtest, mu_train, ...){
  # Step 1: Fit subgroup model #
  mod <- # model call 
  # Step 2: Predictions #
  # Option 1: Create Prediction Function #
  pred.fun <- function(mod, X=NULL){
    Subgrps <- # Predict subgroup assignment
    return( list(Subgrps=Subgrps) )
  }
  # Option 2: Output Subgroups for train/test (here we use pred.fun)
  Subgrps.train = pred.fun(mod, X)
  Subgrps.test = pred.fun(mod, X)
  #Return fit and pred.fun (or just Subgrps.train/Subgrps.test)
  res <- list(mod=mod, pred.fun=pred.fun, Subgrps.train=Subgrps.train,
                  Subgrps.test=Subgrps.test)
  return(res)
}

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
param_template <- function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s,...){
  # Key Outputs: Subgroup specific and overall parameter estimates
  # Overall/Subgroup Specific Estimate ##
  looper = function(s, alpha){
    # Extract parameter estimates #
    return( summ )
  }
   # Across Subgroups #
  S_levels = as.numeric( names(table(Subgrps)) )
  param.dat = lapply(S_levels, looper, alpha_s)
  param.dat = do.call(rbind, param.dat)
  param.dat = data.frame( param.dat )
  ## Overall ##
  param.dat0 = looper(S_levels, alpha_ovrl)
  # Combine and return ##
  param.dat = rbind(param.dat0, param.dat)
  return( param.dat )
}

## ----user_param----------------------------------------------------------

### Robust linear Regression: E(Y|A=1) - E(Y|A=0) ###
param_rlm = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, ...){
  require(MASS)
  indata = data.frame(Y=Y,A=A, X)

  ## Subgroup Specific Estimate ##
  looper = function(s, alpha){
    rlm.mod = tryCatch( rlm(Y ~ A , data=indata[Subgrps %in% s,]),
                       error = function(e) "param error" )
    n.s = dim(indata[Subgrps %in% s,])[1]
    est = summary(rlm.mod)$coefficients[2,1]
    SE = summary(rlm.mod)$coefficients[2,2]
    LCL =  est-qt(1-alpha/2, n.s-1)*SE
    UCL =  est+qt(1-alpha/2, n.s-1)*SE
    pval = 2*pt(-abs(est/SE), df=n.s-1)
    summ <- data.frame(estimand = "E(Y|A=1)-E(Y|A=0)", 
                       Subgrps = ifelse(n.s==dim(X)[1], 0, s),
                       N= n.s, est=est, SE=SE, LCL=LCL, UCL=UCL, pval=pval)
    return( summ )
  }
  # Across Subgroups #
  S_levels = as.numeric( names(table(Subgrps)) )
  param.dat = lapply(S_levels, looper, alpha_s)
  param.dat = do.call(rbind, param.dat)
  param.dat = data.frame( param.dat )
  ## Overall ##
  param.dat0 = looper(S_levels, alpha_ovrl)
  # Combine and return ##
  param.dat = rbind(param.dat0, param.dat)
  return( param.dat )
}


## ----user_prism, warnings=FALSE, message=FALSE---------------------------

res_user1 = PRISM(Y=Y, A=A, X=X, family="gaussian", filter="filter_lasso", 
             ple = "ple_ranger_mtry", submod = "submod_lmtree_pred",
             param="param_rlm")
## variables that remain after filtering ##
res_user1$filter.vars
## Subgroup model: lmtree searching for predictive only ##
plot(res_user1)
## Parameter estimates/inference
res_user1$param.dat
## Waterfall plot of individual treatment effects
plot(res_user1, type="PLE:waterfall")


