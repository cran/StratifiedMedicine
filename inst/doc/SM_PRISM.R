## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=10, fig.height=8.5
)

## ----sim_ctns, warning=FALSE, message=FALSE-----------------------------------
library(ggplot2)
library(dplyr)
library(partykit)
library(StratifiedMedicine)
library(survival)
dat_ctns = generate_subgrp_data(family="gaussian")
Y = dat_ctns$Y
X = dat_ctns$X # 50 covariates, 46 are noise variables, X1 and X2 are truly predictive
A = dat_ctns$A # binary treatment, 1:1 randomized 
length(Y)
table(A)
dim(X)

## ----filter_glmnet, warning=FALSE, message=FALSE------------------------------
res_f <- filter_train(Y, A, X, filter="glmnet")
res_f$filter.vars
plot_importance(res_f)

## ----filter_glmnet2, warning=FALSE, message=FALSE, include=FALSE--------------
res_f2 <- filter_train(Y, A, X, filter="glmnet", hyper=list(interaction=T))
res_f2$filter.vars
plot_importance(res_f2)

## ----ple_train, warning=FALSE, message=FALSE----------------------------------
res_p1 <- ple_train(Y, A, X, ple="ranger", meta="X-learner")
summary(res_p1$mu_train)
plot_ple(res_p1, target = "diff_1_0") + 
  ggtitle("Waterfall Plot: E(Y|A=1)-E(Y|A=0)") + ylab("E(Y|A=1)-E(Y|A=0)")
plot_dependence(res_p1, X=X, vars="X1") + ylab("E(Y|A=1)-E(Y|A=0)")

## ----ple_train2, warning=FALSE, message=FALSE---------------------------------
res_p2 <- ple_train(Y, A, X, ple="ranger", meta="T-learner", hyper=list(mtry=5))
summary(res_p2$mu_train)
plot_dependence(res_p2, X=X, vars=c("X1", "X2")) + 
  ggtitle("Heat Map (By X1,X2): E(Y|A=1)-E(Y|A=0)")
  ylab("E(Y|A=1)-E(Y|A=0)")

## ----submod_train1, warning=FALSE, message=FALSE------------------------------
res_s1 <- submod_train(Y, A, X, submod="lmtree")
table(res_s1$Subgrps.train)
plot(res_s1$fit$mod)

## ----submod_train2, warning=FALSE, message=FALSE------------------------------
res_s2 <- submod_train(Y, A, X,  mu_train=res_p2$mu_train, 
                       submod="otr", hyper=list(thres=">0.10"))
plot(res_s2$fit$mod)

## ----param1, warning=FALSE, message=FALSE-------------------------------------
param.dat1 <- param_est(Y, A, X, Subgrps = res_s1$Subgrps.train, param="lm")
param.dat1 %>% filter(estimand=="mu_1-mu_0")

## ----param2, warning=FALSE, message=FALSE-------------------------------------
param.dat2 <- param_est(Y, A, X, Subgrps = res_s1$Subgrps.train, 
                        mu_hat = res_p1$mu_train, param="dr")
param.dat2 %>% filter(estimand=="mu_1-mu_0")

## ----table_steps, echo=FALSE--------------------------------------------------
library(knitr)
summ.table = data.frame( `Step` = c("estimand(s)", "filter", "ple", "submod", "param"),
                        `gaussian` = c("E(Y|A=0)<br>E(Y|A=1)<br>E(Y|A=1)-E(Y|A=0)",
                                       "Elastic Net<br>(glmnet)", 
                                       "X-learner: Random Forest<br>(ranger)",
                                       "MOB(OLS)<br>(lmtree)", 
                                       "Double Robust<br>(dr)"),
                        `binomial` = c("E(Y|A=0)<br>E(Y|A=1)<br>E(Y|A=1)-E(Y|A=0)",
                                       "Elastic Net<br>(glmnet)", 
                                       "X-learner: Random Forest<br>(ranger)",
                                       "MOB(GLM)<br>(glmtree)", 
                                       "Doubly Robust<br>(dr)"),    
                        `survival` = c("HR(A=1 vs A=0)",
                                       "Elastic Net<br>(glmnet)", 
                                       "T-learner: Random Forest<br>(ranger)",
                                       "MOB(weibull)<br>(mob_weib)", 
                                       "Hazard Ratios<br>(cox)") )                        
                      
kable( summ.table, caption = "Default PRISM Configurations (With Treatment)", full_width=T)

summ.table = data.frame(`Step` = c("estimand(s)", "filter", "ple", "submod", "param"),
                        `gaussian` = c("E(Y)",
                                       "Elastic Net<br>(glmnet)", 
                                       "Random Forest<br>(ranger)",
                                       "Conditional Inference Trees<br>(ctree)",
                                       "OLS<br>(lm)"),
                        `binomial` = c("Prob(Y)",
                                       "Elastic Net<br>(glmnet)", 
                                       "Random Forest<br>(ranger)",
                                       "Conditional Inference Trees<br>(ctree)", 
                                       "OLS<br>(lm)"),    
                        `survival` = c("RMST", "Elastic Net<br>(glmnet)", 
                                       "Random Forest<br>(ranger)",
                                       "Conditional Inference Trees<br>(ctree)",
                                       "RMST<br>(rmst)"))                        
                      
kable( summ.table, caption = "Default PRISM Configurations (Without Treatment, A=NULL)", full_width=T)

## ----default_ctns, warning=FALSE----------------------------------------------
# PRISM Default: filter_glmnet, ranger, lmtree, dr #
res0 = PRISM(Y=Y, A=A, X=X)
summary(res0)
plot(res0) # same as plot(res0, type="tree")

## ----default_ctns_prog, warning=FALSE-----------------------------------------
# PRISM Default: filter_glmnet, ranger, ctree, param_lm #
res_prog = PRISM(Y=Y, X=X)
# res_prog = PRISM(Y=Y, A=NULL, X=X) #also works
summary(res_prog)

## ----default_ctns_filter, include=FALSE---------------------------------------
# elastic net model: loss by lambda #
plot(res0$filter.mod)
## Variables that remain after filtering ##
res0$filter.vars
# All predictive variables (X1,X2) and prognostic variables (X3,X5, X7) remains.
plot_importance(res0)

## ----default_ctns_ple---------------------------------------------------------
summary(res0$mu_train)
plot_ple(res0)
plot_dependence(res0, vars=c("X2"))

## ----default_ctns_submod------------------------------------------------------
plot(res0$submod.fit$mod, terminal_panel = NULL)
table(res0$out.train$Subgrps)
table(res0$out.test$Subgrps)

## ----default_ctns2------------------------------------------------------------
## Overall/subgroup specific parameter estimates/inference
res0$param.dat

## ----default_hyper------------------------------------------------------------
# PRISM Default: glmnet, ranger, lmtree, dr #
# Change hyper-parameters #
res_new_hyper = PRISM(Y=Y, A=A, X=X, filter.hyper = list(lambda="lambda.1se"),
                      ple.hyper = list(min.node.pct=0.05), 
                      submod.hyper = list(minsize=200), verbose=FALSE)
summary(res_new_hyper)

## ----default_binary-----------------------------------------------------------
dat_bin = generate_subgrp_data(family="binomial", seed = 5558)
Y = dat_bin$Y
X = dat_bin$X # 50 covariates, 46 are noise variables, X1 and X2 are truly predictive
A = dat_bin$A # binary treatment, 1:1 randomized 

res0 = PRISM(Y=Y, A=A, X=X)
summary(res0)
plot(res0)

## ----default_surv-------------------------------------------------------------
# Load TH.data (no treatment; generate treatment randomly to simulate null effect) ##
data("GBSG2", package = "TH.data")
surv.dat = GBSG2
# Design Matrices ###
Y = with(surv.dat, Surv(time, cens))
X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
set.seed(6345)
A = rbinom(n = dim(X)[1], size=1, prob=0.5)

# Default: glmnet ==> ranger (estimates patient-level RMST(1 vs 0) ==> mob_weib (MOB with Weibull) ==> cox (Cox regression)
res_weib = PRISM(Y=Y, A=A, X=X)
plot(res_weib, type="PLE:waterfall")
plot(res_weib)

## ----default_boot, warning=FALSE, message=FALSE-------------------------------
res_boot = PRISM(Y=Y, A=A, X=X, resample = "Bootstrap", R=50, ple="None")
summary(res_boot)
# Plot of distributions #
plot(res_boot, type="resample", estimand = "HR(A=1 vs A=0)")+geom_vline(xintercept = 1)

