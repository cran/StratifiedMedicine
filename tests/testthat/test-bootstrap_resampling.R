test_that("Bootstrap resampling produces expected resamples per subgroup (no pooling)", {

  skip_on_cran()

  library(partykit)

  dat_ctns <- generate_subgrp_data(family = "gaussian")
  Y <- dat_ctns$Y
  X <- dat_ctns$X
  A <- dat_ctns$A

  R <- 30
  res <- PRISM(Y = Y, A = A, X = X, ple = "None", param = "lm",
               resample = "Bootstrap", R = R)

  # resamp_dist should exist
  expect_false(is.null(res$resamp_dist))

  rdist <- res$resamp_dist
  subgrps <- unique(res$out.train$Subgrps)
  estimands <- unique(rdist$estimand)

  # For each subgroup + overall, check we get close to R resamples
  for (s in c("ovrl", subgrps)) {
    for (e in estimands) {
      n_resamples <- sum(rdist$Subgrps == s & rdist$estimand == e, na.rm = TRUE)
      expect_gte(n_resamples, floor(R * 0.90),
                 label = paste("Resamples for Subgrps =", s, ", estimand =", e))
    }
  }

  # param.dat should contain treatment estimates (bootstrap-adjusted)
  expect_false(is.null(res$param.dat))
  expect_true(nrow(res$param.dat) > 0)
})

test_that("Bootstrap resampling with pooling produces expected resamples", {

  skip_on_cran()

  library(partykit)

  dat_ctns <- generate_subgrp_data(family = "gaussian")
  Y <- dat_ctns$Y
  X <- dat_ctns$X
  A <- dat_ctns$A

  R <- 30
  res <- PRISM(Y = Y, A = A, X = X, ple = "None", param = "lm",
               pool = "trteff", resample = "Bootstrap", R = R)

  # resamp_dist should exist
  expect_false(is.null(res$resamp_dist))

  rdist <- res$resamp_dist
  estimands <- unique(rdist$estimand)

  # For overall, check we get close to R resamples
  for (e in estimands) {
    n_ovrl <- sum(rdist$Subgrps == "ovrl" & rdist$estimand == e, na.rm = TRUE)
    expect_gte(n_ovrl, floor(R * 0.90),
               label = paste("Overall resamples for estimand =", e))
  }

  # With pooling, trt_assign should exist
  expect_false(is.null(res$trt_assign))

  # param.dat should contain treatment estimates (bootstrap-adjusted)
  expect_false(is.null(res$param.dat))
  expect_true(nrow(res$param.dat) > 0)
})

test_that("Bootstrap resampling works for survival data", {

  skip_on_cran()

  library(survival)
  require(TH.data)

  data("GBSG2", package = "TH.data")
  surv.dat <- GBSG2
  Y <- with(surv.dat, Surv(time, cens))
  X <- surv.dat[, !(colnames(surv.dat) %in% c("time", "cens"))]
  set.seed(513)
  A <- rbinom(n = nrow(X), size = 1, prob = 0.5)

  R <- 20
  res <- PRISM(Y = Y, A = A, X = X, resample = "Bootstrap", R = R,
               ple = "None")

  expect_false(is.null(res$resamp_dist))

  rdist <- res$resamp_dist
  estimands <- unique(rdist$estimand)

  # Overall should have close to R resamples per estimand
  for (e in estimands) {
    n_ovrl <- sum(rdist$Subgrps == "ovrl" & rdist$estimand == e, na.rm = TRUE)
    expect_gte(n_ovrl, floor(R * 0.90),
               label = paste("Survival overall resamples for estimand =", e))
  }
})
