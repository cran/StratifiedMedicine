test_that("Test different resampling approaches (continuous)", {
  
  skip_on_cran()
  
  ## Continuous ##
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  
  # ### Bootstrap ###
  res0 <- PRISM(Y, A, X, ple="None", param="lm", resample="Bootstrap", R=50)
  test0 = ifelse( is.null(res0$resamp_dist), 0, 1)
  
  ### Cross-Validation ###
  res1 <- PRISM(Y, A, X, ple="None", param="lm", resample="CV")
  test1 = ifelse( is.null(res1$resamp_dist), 0, 1)

  tests_ctns = ifelse( sum(test0,test1)==2, 1, 0)

  ### Output Test Results ###
  expect_equal(tests_ctns, 1L)
})

test_that("Test different resampling approaches (survival)", {
  skip_on_cran()
  
  ### Survival Tests ###
  library(survival)
  require(TH.data); require(coin)
  data("GBSG2", package = "TH.data")
  surv.dat = GBSG2
  # Design Matrices ###
  Y = with(surv.dat, Surv(time, cens))
  X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
  set.seed(513)
  A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
  
  ### Bootstrap ###
  res0 <- PRISM(Y, A, X, resample="Bootstrap", R=50, ple="None")
  test0 = ifelse( is.null(res0$resamp_dist), 0, 1)
  
  ### Cross-Validation ###
  res1 <- PRISM(Y, A, X, resample="CV", ple="None")
  test1 = ifelse( is.null(res1$resamp_dist), 0, 1)
  
  tests_surv = ifelse( sum(test0,test1)==2, 1, 0)
  
  ### Output Test Results ###
  expect_equal(tests_surv, 1L)
})