## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=8, fig.height=6 
)

## ----sim_ctns------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(partykit)
library(StratifiedMedicine)
dat_ctns = generate_subgrp_data(family="gaussian", seed=65532)
Y = dat_ctns$Y
X = dat_ctns$X # 50 covariates, 46 are noise variables, X1 and X2 are truly predictive
A = dat_ctns$A # binary treatment, 1:1 randomized 
length(Y)
table(A)
dim(X)

## ----default_ctns--------------------------------------------------------
# PRISM Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
res0 = PRISM(Y=Y, A=A, X=X)
## This is the same as running ##
# res1 = PRISM(Y=Y, A=A, X=X, family="gaussian", filter="filter_glmnet", 
#              ple = "ple_ranger", submod = "submod_lmtree", param="param_ple")

## ----default_ctns_filter-------------------------------------------------
# elastic net model: loss by lambda #
plot(res0$filter.mod)
## Variables that remain after filtering ##
res0$filter.vars
# All predictive variables (X1,X2) and prognostic variables (X3,X5, X7) remains.

## ----default_ctns_ple----------------------------------------------------
prob.PLE = mean(I(res0$mu_train$PLE>0))
# Density Plot #
plot(res0, type="PLE:density")+geom_vline(xintercept = 0) +
     geom_text(x=1.5, y=0.4, label=paste("Prob(PLE>0)=", prob.PLE, sep=""))
# Waterfall Plot #
plot(res0, type="PLE:waterfall")+geom_vline(xintercept = 0) + 
  geom_text(x=200, y=1, label=paste("Prob(PLE>0)=", prob.PLE, sep=""))

## ----default_ctns_submod-------------------------------------------------
plot(res0$submod.fit$mod, terminal_panel = NULL)
table(res0$out.train$Subgrps)
table(res0$out.test$Subgrps)

## ----default_ctns2-------------------------------------------------------
## Overall/subgroup specific parameter estimates/inference
res0$param.dat
## Forest plot: Overall/subgroup specific parameter estimates (CIs)
plot(res0, type="forest")

## ----modify_submod_plot--------------------------------------------------
param.dat = res0$param.dat[res0$param.dat$Subgrps>0,]
param.dat$est = with(param.dat, sprintf("%.2f",round(est,2)))
param.dat$CI = with(param.dat, paste("[", sprintf("%.2f",round(param.dat$LCL,2)),",",
                                     sprintf("%.2f",round(param.dat$UCL,2)),"]",sep=""))
mu_1 = aggregate(res0$mu_train$mu1~res0$out.train$Subgrps, FUN="mean")
colnames(mu_1) = c("Subgrps", "est.A1")
mu_0 = aggregate(res0$mu_train$mu0~res0$out.train$Subgrps, FUN="mean")
colnames(mu_0) = c("Subgrps", "est.A0")
param.dat = left_join(param.dat, mu_1, by="Subgrps")
param.dat = left_join(param.dat, mu_0, by="Subgrps")
param.dat$est.A1 = with(param.dat, sprintf("%.2f",round(est.A1,2)))
param.dat$est.A0 = with(param.dat, sprintf("%.2f",round(est.A0,2)))
smod = res0$submod.fit$mod
smod_node <- as.list(smod$node)
for(i in 1:nrow(param.dat)){
   smod_node[[param.dat[i,1]]]$info$est <- param.dat$est[i]
   smod_node[[param.dat[i,1]]]$info$CI <-  param.dat$CI[i]
   smod_node[[param.dat[i,1]]]$info$est.A1 <- param.dat$est.A1[i]
   smod_node[[param.dat[i,1]]]$info$est.A0 <- param.dat$est.A0[i]
}
smod$node <- as.partynode(smod_node)
plot(smod, terminal_panel = node_terminal, tp_args = list(
  FUN = function(node) c( paste("n =", node$nobs),
                          paste("E(Y|A=0):", node$est.A0),
                          paste("E(Y|A=1):", node$est.A1),
                          paste("Diff:",node$est),
                          node$CI) ) )


## ----heat_maps-----------------------------------------------------------
grid.data = expand.grid(X1 = seq(min(X$X1), max(X$X1), by=0.30),
                    X2 = seq(min(X$X2), max(X$X2), by=0.30))
plot(res0, type="heatmap", grid.data = grid.data)


## ----default_hyper-------------------------------------------------------
# PRISM Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
# Change hyper-parameters #
res_new_hyper = PRISM(Y=Y, A=A, X=X, filter.hyper = list(lambda="lambda.1se"),
                      ple.hyper = list(min.node.pct=0.05), 
                      submod.hyper = list(minsize=200))
plot(res_new_hyper$submod.fit$mod) # Plot subgroup model results
plot(res_new_hyper) # Forest plot 

## ----default_boot--------------------------------------------------------
library(ggplot2)
library(dplyr)
res_boot = PRISM(Y=Y, A=A, X=X, resample = "Bootstrap", R=50, verbose=FALSE)
# # Plot of distributions and P(est>0) #
plot(res_boot, type="resample")+geom_vline(xintercept = 0)
aggregate(I(est>0)~Subgrps, data=res_boot$resamp.dist, FUN="mean")


## ----default_surv--------------------------------------------------------
library(survival)
library(ggplot2)
# Load TH.data (no treatment; generate treatment randomly to simulate null effect) ##
data("GBSG2", package = "TH.data")
surv.dat = GBSG2
# Design Matrices ###
Y = with(surv.dat, Surv(time, cens))
X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
A = rbinom( n = dim(X)[1], size=1, prob=0.5  )

# Default: filter_glmnet ==> submod_weibull (MOB with Weibull) ==> param_cox (Cox regression)
res_weibull1 = PRISM(Y=Y, A=A, X=X)
plot(res_weibull1, type="PLE:waterfall")
plot(res_weibull1$submod.fit$mod)
plot(res_weibull1)+ylab("HR [95% CI]")

# PRISM: filter_glmnet ==> submod_ctree ==> param_cox (Cox regression) #
res_ctree1 = PRISM(Y=Y, A=A, X=X, ple=NULL, submod = "submod_ctree")
plot(res_ctree1$submod.fit$mod)
plot(res_ctree1)+ylab("HR [95% CI]")


## ----sim_ctns2-----------------------------------------------------------
library(StratifiedMedicine)
dat_ctns = generate_subgrp_data(family="gaussian")
Y = dat_ctns$Y
X = dat_ctns$X # 50 covariates, 46 are noise variables, X1 and X2 are truly predictive
A = dat_ctns$A # binary treatment, 1:1 randomized 

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
    mods = list(mod0=mod0, mod1=mod1)
    ## Predictions: Train/Test ##
    mu_train = data.frame( mu1 =  predict(mod1, data = X)$predictions,
                             mu0 = predict(mod0, data = X)$predictions)
    mu_train$PLE = with(mu_train, mu1 - mu0 )

    mu_test = data.frame( mu1 =  predict(mod1, data = Xtest)$predictions,
                            mu0 = predict(mod0, data = Xtest)$predictions)
    mu_test$PLE = with(mu_test, mu1 - mu0 )
    res = list(mods=mods, mu_train=mu_train, mu_test=mu_test)
    class(res) = "ple_ranger_mtry"
    return( list(mods=mods, mu_train=mu_train, mu_test=mu_test))
}

## ----user_submod---------------------------------------------------------
submod_lmtree_pred = function(Y, A, X, Xtest, mu_train, ...){
  require(partykit)
  ## Fit Model ##
  mod <- lmtree(Y~A | ., data = X, parm=2) ##parm=2 focuses on treatment interaction #
  ##  Predict Subgroups for Train/Test ##
  Subgrps.train = as.numeric( predict(mod, type="node") )
  Subgrps.test = as.numeric( predict(mod, type="node", newdata = Xtest) )
  ## Predict E(Y|X=x, A=1)-E(Y|X=x,A=0) ##
  pred.train = predict( mod, newdata = data.frame(A=1, X) ) -
    predict( mod, newdata = data.frame(A=0, X) )
  pred.test =  predict( mod, newdata = data.frame(A=1, Xtest) ) -
    predict( mod, newdata = data.frame(A=0, Xtest) )
  ## Return Results ##
  return(  list(mod=mod, Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test,
                pred.train=pred.train, pred.test=pred.test) )
}

## ----user_param----------------------------------------------------------

### Robust linear Regression: E(Y|A=1) - E(Y|A=0) ###
param_rlm = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, ...){
  require(MASS)
  indata = data.frame(Y=Y,A=A, X)
  mod.ovrl = rlm(Y ~ A , data=indata)
  param.dat0 = data.frame( Subgrps=0, N = dim(indata)[1],
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
  param.dat = data.frame( S = S_levels, N=S_N, param.dat)
  colnames(param.dat) = c("Subgrps", "N", "est", "SE", "LCL", "UCL", "pval")
  param.dat = rbind( param.dat0, param.dat)
  return( param.dat )
}


## ----user_prism----------------------------------------------------------

res_user1 = PRISM(Y=Y, A=A, X=X, family="gaussian", filter="filter_lasso", 
             ple = "ple_ranger_mtry", submod = "submod_lmtree_pred",
             param="param_rlm")
## variables that remain after filtering ##
res_user1$filter.vars
## Subgroup model: lmtree searching for predictive only ##
plot(res_user1$submod.fit$mod)
## Parameter estimates/inference
res_user1$param.dat
## Forest Plot (95% CI) ##
plot(res_user1)


