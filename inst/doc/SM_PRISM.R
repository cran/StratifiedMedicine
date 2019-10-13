## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=11.5, fig.height=8.5
)

## ----table_steps, echo=FALSE---------------------------------------------
library(knitr)
summ.table = data.frame( `Step` = c("estimand(s)", "filter", "ple", "submod", "param"),
                        `gaussian` = c("E(Y|A=0)<br>E(Y|A=1)<br>E(Y|A=1)-E(Y|A=0)",
                                       "filter_glmnet", "ple_ranger",
                                       "submod_lmtree", "param_ple"),
                        `binomial` = c("E(Y|A=0)<br>E(Y|A=1)<br>E(Y|A=1)-E(Y|A=0)",
                                       "filter_glmnet", "ple_ranger",
                                       "submod_lmtree", "param_ple"),    
                        `survival` = c("HR(A=1 vs A=0)", "filter_glmnet", "ple_ranger",
                                       "submod_weibull", "param_ple") )                        
                      
kable( summ.table, caption = "Default PRISM Configurations (With Treatment)", full_width=T)

summ.table = data.frame( `Step` = c("estimand(s)", "filter", "ple", "submod", "param"),
                        `gaussian` = c("E(Y)",
                                       "filter_glmnet", "ple_ranger",
                                       "submod_ctree", "param_lm"),
                        `binomial` = c("Prob(Y)",
                                       "filter_glmnet", "ple_ranger",
                                       "submod_ctree", "param_lm"),    
                        `survival` = c("RMST", "filter_glmnet", "ple_ranger",
                                       "submod_ctree", "param_rmst") )                        
                      
kable( summ.table, caption = "Default PRISM Configurations (Without Treatment, A=NULL)", full_width=T)

## ----sim_ctns, warning=FALSE, message=FALSE------------------------------
library(ggplot2)
library(dplyr)
library(partykit)
library(StratifiedMedicine)
dat_ctns = generate_subgrp_data(family="gaussian")
Y = dat_ctns$Y
X = dat_ctns$X # 50 covariates, 46 are noise variables, X1 and X2 are truly predictive
A = dat_ctns$A # binary treatment, 1:1 randomized 
length(Y)
table(A)
dim(X)

## ----default_ctns, warning=FALSE-----------------------------------------
# PRISM Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
res0 = PRISM(Y=Y, A=A, X=X)
summary(res0)
plot(res0) # same as plot(res0, type="submod")
## This is the same as running ##
# res1 = PRISM(Y=Y, A=A, X=X, family="gaussian", filter="filter_glmnet", 
#              ple = "ple_ranger", submod = "submod_lmtree", param="param_ple")

## ----default_ctns_prog, warning=FALSE------------------------------------
# PRISM Default: filter_glmnet, ple_ranger, submod_ctree, param_lm #
res_prog = PRISM(Y=Y, X=X)
summary(res_prog)
plot(res_prog)

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
     geom_text(x=1, y=0.4, label=paste("Prob(PLE>0)=", prob.PLE, sep=""))
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
plot(res0, type="submod")
plot(res0, type="forest")

## ----heat_maps-----------------------------------------------------------
grid.data = expand.grid(X1 = seq(min(X$X1), max(X$X1), by=0.30),
                    X2 = seq(min(X$X2), max(X$X2), by=0.30))
plot(res0, type="heatmap", grid.data = grid.data)


## ----default_hyper-------------------------------------------------------
# PRISM Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
# Change hyper-parameters #
res_new_hyper = PRISM(Y=Y, A=A, X=X, filter.hyper = list(lambda="lambda.1se"),
                      ple.hyper = list(min.node.pct=0.05), 
                      submod.hyper = list(minsize=200), verbose=FALSE)
plot(res_new_hyper)

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

# Default: filter_glmnet ==> ple_ranger (estimates patient-level RMST(1 vs 0) ==> submod_weibull (MOB with Weibull) ==> param_cox (Cox regression)
res_weibull1 = PRISM(Y=Y, A=A, X=X)
plot(res_weibull1, type="PLE:waterfall")
plot(res_weibull1)

# PRISM: filter_glmnet ==> submod_ctree ==> param_cox (Cox regression) #
res_ctree1 = PRISM(Y=Y, A=A, X=X, submod = "submod_ctree")
plot(res_ctree1)


## ----default_boot, warning=FALSE, message=FALSE--------------------------
library(ggplot2)
library(dplyr)
res_boot = PRISM(Y=Y, A=A, X=X, resample = "Bootstrap", R=50, ple = "None")
# Plot of distributions #
plot(res_boot, type="resample", estimand = "HR(A=1 vs A=0)")+geom_vline(xintercept = 1)

