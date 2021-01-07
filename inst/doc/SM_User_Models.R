## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=9.5, fig.height=8.5 
)

## ----sim_ctns0----------------------------------------------------------------
library(StratifiedMedicine)
dat_ctns = generate_subgrp_data(family="gaussian")
Y = dat_ctns$Y
X = dat_ctns$X # 50 covariates, 46 are noise variables, X1 and X2 are truly predictive
A = dat_ctns$A # binary treatment, 1:1 randomized 

## ----table_steps, echo=FALSE--------------------------------------------------
library(knitr)
summ.table = data.frame( `Model` = c("filter", "ple", "submod", "param"),
                        `Required Outputs` = c("filter.vars", "list(mod,pred.fun)",
                                               "list(mod,pred.fun)", 
                                               "param.dat"),
                        `Description` = c("Variables that pass filter",
                                          "Model fit(s) and prediction function",
                                          "Model fit(s) and prediction function",
                                          "Parameter Estimates") )
kable( summ.table, caption = "Key Outputs by Model" )

## ----user_filter_template-----------------------------------------------------
filter_template = function(Y, A, X, ...){
  # Step 1: Fit Filter Model #
  mod <- # model call 
  # Step 2: Extract variables that pass the filter #
  filter.vars <- # depends on mod fit
  # Return model fit and filtered variables #
  res = list(mod=mod, filter.vars=filter.vars)
  return( res )
}

## ----user_filter--------------------------------------------------------------
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

## ----user_ple_template--------------------------------------------------------
ple_template <- function(Y, A, X, ...){
  # Step 1: Fit PLE Model #
  # for example: Estimate E(Y|A=1,X), E(Y|A=0,X), E(Y|A=1,X)-E(Y|A=0,X)
  mod <- # ple model call 
  # mod = list(mod0=mod0, mod1=mod1) # If multiple fitted models, combine into list
  # Step 2: Predictions
  # Option 1: Create a Prediction Function #
  pred.fun <- function(mod, X, ...){
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

## ----user_ple-----------------------------------------------------------------
ple_ranger_mtry = function(Y, X, mtry=5, ...){
   require(ranger)
    train =  data.frame(Y=Y, X)
    mod <- ranger(Y ~ ., data = train, seed=1, mtry = mtry)
    mod = list(mod=mod)
    pred.fun <- function(mod, X, ...){
      mu_hat <- predict(mod$mod, X)$predictions
      mu_hat <- data.frame(mu_hat)
      return(mu_hat)
    }
    res = list(mod=mod, pred.fun=pred.fun)
    return(res)
}

## ----user_submod_template-----------------------------------------------------
submod_template <- function(Y, A, X, Xtest, mu_train, ...){
  # Step 1: Fit subgroup model #
  mod <- # model call 
  # Step 2: Predictions #
  # Option 1: Create Prediction Function #
  pred.fun <- function(mod, X=NULL, ...){
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

## ----user_submod--------------------------------------------------------------
submod_lmtree_pred = function(Y, A, X, mu_train, ...){
  require(partykit)
  ## Fit Model ##
  mod <- lmtree(Y~A | ., data = X, parm=2) ##parm=2 focuses on treatment interaction #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
     Subgrps <- NULL
     Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
     return( list(Subgrps=Subgrps) )
  }
  ## Return Results ##
  return(list(mod=mod, pred.fun=pred.fun))
}

## ----user_param_template------------------------------------------------------
param_template <- function(Y, A, X, mu_hat, alpha,...){
  # Key Outputs: Subgroup specific and overall parameter estimates
  mod <- # Call parameter model #
  # Extract estimates/variability and combine #
  param.dat <- data.frame(n=n, estimand="mu_1", 
                          est=est, SE=SE, LCL=LCL, UCL=UCL, pval=pval)
  return(param.dat)
}

## ----user_param---------------------------------------------------------------
### Robust linear Regression: E(Y|A=1) - E(Y|A=0) ###
param_rlm = function(Y, A, alpha, ...){
  require(MASS)
  indata = data.frame(Y=Y,A=A)
  rlm.mod = tryCatch( rlm(Y ~ A , data=indata),
                       error = function(e) "param error" )
  n = dim(indata)[1]
  est = summary(rlm.mod)$coefficients[2,1]
  SE = summary(rlm.mod)$coefficients[2,2]
  LCL =  est-qt(1-alpha/2, n-1)*SE
  UCL =  est+qt(1-alpha/2, n-1)*SE
  pval = 2*pt(-abs(est/SE), df=n-1)
  param.dat <- data.frame(N= n, estimand = "mu_1-mu_0",
                     est=est, SE=SE, LCL=LCL, UCL=UCL, pval=pval)
  return(param.dat)
}


## ----user_SM_final, warnings=FALSE, message=FALSE-----------------------------
step1 <- filter_train(Y, A, X, filter="filter_lasso")
X.star <- X[,colnames(X) %in% step1$filter.vars]
step2 <- ple_train(Y, A, X.star, ple = "ple_ranger_mtry")
plot_ple(step2)
step3 <- submod_train(Y, A, X.star, submod = "submod_lmtree_pred")
plot(step3$fit$mod)
step4 <- param_est(Y, A, X.star, Subgrps=step3$Subgrps.train, param="param_rlm")
step4
# PRISM #
res_user1 = PRISM(Y=Y, A=A, X=X, family="gaussian", filter="filter_lasso", 
             ple = "ple_ranger_mtry", submod = "submod_lmtree_pred",
             param="param_rlm")
plot_ple(res_user1)
# variables that remain after filtering #
res_user1$filter.vars
# Parameter estimates/inference
res_user1$param.dat

