% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ple_train.R
\name{ple_train}
\alias{ple_train}
\title{Patient-level Estimates: Train Model}
\usage{
ple_train(
  Y,
  A,
  X,
  Xtest = NULL,
  family = "gaussian",
  propensity = FALSE,
  ple = "ranger",
  meta = ifelse(family == "survival", "T-learner", "X-learner"),
  hyper = NULL,
  tau = NULL,
  ...
)
}
\arguments{
\item{Y}{The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )}

\item{A}{Treatment variable. (Default supports binary treatment, either numeric or 
factor). "ple_train" accomodates >2 along with binary treatments.}

\item{X}{Covariate space.}

\item{Xtest}{Test set. Default is NULL (no test predictions). Variable types should match X.}

\item{family}{Outcome type. Options include "gaussion" (default), "binomial", and "survival".}

\item{propensity}{Propensity score estimation, P(A=a|X). Default=FALSE which 
use the marginal estimates, P(A=a) (applicable for RCT data). If TRUE, will 
use the "ple" base learner to estimate P(A=a|X).}

\item{ple}{Base-learner used to estimate patient-level equantities, such as the 
conditional average treatment effect (CATE), E(Y|A=1,X)-E(Y|A=0, X) = CATE(X). 
Default is random based based through "ranger". "None" uses no ple. See below for 
details on estimating the treatment contrasts.}

\item{meta}{Using the ple model as a base learner, meta-learners can be used for 
estimating patient-level treatment differences. Options include "T-learner" (treatment
specific models), "S-learner" (single model), and "X-learner". For family="gaussian" &
"binomial", the default is "X-learner", which uses a two-stage regression 
approach (See Kunzel et al 2019). For "survival", the default is "T-learner". "X-learner" 
is currently not supported for survival outcomes.}

\item{hyper}{Hyper-parameters for the ple model (must be list). Default is NULL.}

\item{tau}{Maximum follow-up time for RMST based estimates (family="survival"). 
Default=NULL, which takes min(max(time[a])), for a=1,..,A.}

\item{...}{Any additional parameters, not currently passed through.}
}
\value{
Trained ple models and patient-level estimates for train/test sets. 
 \itemize{
  \item mod - trained model(s)
  \item mu_train - Patient-level estimates (training set)
  \item mu_test - Patient-level estimates (test set)
}
}
\description{
Wrapper function to train a patient-level estimate (ple) model. Used directly in PRISM
and can be used to directly fit a ple model by name.
}
\details{
ple_train uses base-learners along with a meta-learner to obtain patient-level 
estimates under different treatment exposures (see Kunzel et al).  For family="gaussian" 
or "binomial", output estimates of \eqn{\mu(a,x)=E(Y|x,a)} and treatment differences 
(average treatment effect or risk difference). For survival, either logHR based estimates
or RMST based estimates can be obtained. Current base-learner ("ple") options include:


1. \strong{linear}: Uses either linear regression (family="gaussian"), 
logistic regression (family="binomial"), or cox regression (family="survival"). 
No hyper-parameters.

2. \strong{ranger}: Uses random forest ("ranger" R package). The default hyper-parameters are: 
hyper = list(mtry=NULL, min.node.pct=0.10)

where mtry is number of randomly selected variables (default=NULL; sqrt(dim(X)))
and min.node.pct is the minimum node size as a function of the total data size 
(ex: min.node.pct=10\% requires at least 10% of the data at each node)

3. \strong{glmnet}: Uses elastic net ("glmnet" R package). The default hyper-parameters are: 
hyper = list(lambda="lambda.min")

where lambda controls the penalty parameter for predictions. lambda="lambda.1se"
will likely result in a less complex model. 

4. \strong{bart}:  Uses bayesian additive regression trees (Chipman et al 2010; 
BART R package). Default hyper-parameters are:

hyper = list(sparse=FALSE)

where sparse controls whether to perform variable selection based on a sparse 
Dirichlet prior rather than simply uniform.
}
\examples{
\donttest{
library(StratifiedMedicine)
## Continuous ##
dat_ctns = generate_subgrp_data(family="gaussian")
Y = dat_ctns$Y
X = dat_ctns$X
A = dat_ctns$A


# X-Learner (With ranger based learners)
mod1 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="ranger", method="X-learner")
summary(mod1$mu_train)

# T-Learner (Treatment specific)
mod2 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="ranger", method="T-learner")
summary(mod2$mu_train)


mod3 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="bart", method="X-learner")
summary(mod3$mu_train)
}


}
\references{
Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of 
random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. 
\doi{10.18637/jss.v077.i01}

Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for
 Generalized Linear Models via Coordinate Descent,
 \url{https://web.stanford.edu/~hastie/Papers/glmnet.pdf} Journal of Statistical 
 Software, Vol. 33(1), 1-22 Feb 2010 Vol. 33(1), 1-22 Feb 2010.

Chipman, H., George, E., and McCulloch R. (2010) Bayesian Additive 
Regression Trees. The Annals of Applied Statistics, 4,1, 266-298

Kunzel S, Sekhon JS, Bickel PJ, Yu B. Meta-learners for Estimating
Hetergeneous Treatment Effects using Machine Learning. 2019.
}
\seealso{
\code{\link{PRISM}}
}
