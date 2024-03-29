% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param_est.R
\name{param_est}
\alias{param_est}
\title{Parameter Estimation: Across Subgroups}
\usage{
param_est(
  Y,
  A,
  X,
  param,
  mu_hat = NULL,
  Subgrps,
  alpha_ovrl = 0.05,
  alpha_s = 0.05,
  combine = "SS",
  ...
)
}
\arguments{
\item{Y}{The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )}

\item{A}{Treatment variable. (Default supports binary treatment, either numeric or 
factor). "ple_train" accomodates >2 along with binary treatments.}

\item{X}{Covariate space.}

\item{param}{Parameter estimation and inference function. Based on the discovered 
subgroups, estimate parameter estimates and correspond variability metrics. Options
include "lm" (unadjusted linear regression), "dr" (doubly-robust estimator),
"gcomp" (G-computation, average the patient-level estimates), "cox" (cox regression),
and "rmst" (RMST based estimates as in survRMST package). Default for "gaussian",
"binomial" is "dr", while default for "survival" is "cox". Currently only available 
for binary treatments or A=NULL.}

\item{mu_hat}{Patient-level estimates (see \code{ple_train})}

\item{Subgrps}{Identified subgroups. Can be pre-specified, or determined
adaptively (see \code{submod_train}).}

\item{alpha_ovrl}{Two-sided alpha level for overall population}

\item{alpha_s}{Two-sided alpha level at subgroup}

\item{combine}{Given identified subgroups and correspond point-estimates/SEs/sample sizes,
combine="SS" will use sample size weighting for estimates at the overall level. Not 
applicable for param="dr","ple".}

\item{...}{Any additional parameters, not currently passed through.}
}
\value{
Data-set with parameter estimates and corresponding
variability metrics, for overall and subgroups. Subgrps="ovrl" corresponds to the overall
population by default.
 \itemize{
  \item param.dat - Parameter estimates and variability metrics (est, SE,
  LCL/UCL = lower/upper confidence limits, pval = p-value).
  }
}
\description{
For each identified subgroup, obtain point-estimates and variability metrics
(est, SE, CI). fit separate linear regression models. Point-estimates and
variability metrics in the overall population are obtained by aggregating subgroup
specific results (adaptive weighting or sample size weighting).
}
\examples{
library(StratifiedMedicine)

## Continuous ##
dat_ctns = generate_subgrp_data(family="gaussian")
Y = dat_ctns$Y
X = dat_ctns$X
A = dat_ctns$A

## Identify Subgroups: MOB (lmtree) ##
res_lmtree = submod_train(Y, A, X, submod="lmtree")

## Parameter-estimation ##
param.dat = param_est(Y, A, X, param="lm", Subgrps = res_lmtree$Subgrps.train)
param.dat

}
\references{
Funk et al. Doubly Robust Estimation of Causal Effects. 
Am J Epidemiol 2011. 173(7): 761-767.

Andersen, P. and Gill, R. (1982). Cox’s regression model for counting 
processes, a large sample study. Annals of Statistics 10, 1100-1120.

Uno et al. Moving beyond the hazard ratio in quantifying the 
between-group difference in survival analysis. Journal of clinical Oncology 2014, 
32, 2380-2385.
}
\seealso{
\code{\link{param_combine}}
}
