

## Bayesian: Normal prior (at overall for subgroups), normal posterior ##
norm_norm <- function(PRISM.fit, alpha_ovrl, alpha_s,
                      scale = length(PRISM.fit$Subgrps.train), ...) {
  param.dat <- PRISM.fit$param.dat
  Subgrps <- PRISM.fit$Subgrps.train
  get_posterior_ests <- function(param.dat, e) {
    # Set prior to overall #
    param.e <- param.dat[param.dat$estimand==e,]
    prior.mu <- param.e$est[param.e$Subgrps==0]
    prior.var <- scale * ( param.e$SE[param.e$Subgrps==0] )^2
    looper <- function(s){
      ns <- param.e$N[param.e$Subgrps==s]
      obs.mu <- param.e$est[param.e$Subgrps==s]
      obs.var <- ( param.e$SE[param.e$Subgrps==s] )^2
      ### Mean/Var ##
      if (s==0){
        mu_B <- obs.mu
        var_B <- obs.var
        alpha <- alpha_ovrl
      }
      if (s>0){
        var_B <- ( 1/prior.var + 1/(obs.var) )^(-1)
        mu_B <- var_B * ( prior.mu/prior.var + obs.mu / (obs.var)    ) 
        alpha <- alpha_s
      }
      post.prob <- 1-pnorm(0, mean = mu_B, sd = sqrt(var_B) )
      # Bayesian Intervals (q25, q75) #
      LCL <- qnorm( alpha/2, mean = mu_B, sd = sqrt(var_B) )
      UCL <- qnorm( 1-alpha/2, mean = mu_B, sd = sqrt(var_B) )
      summ <- data.frame(Subgrps=s, N=ns, estimand=e,
                        est = mu_B, SE = sqrt(var_B), 
                        LCL=LCL, UCL = UCL, prob.ge0 = post.prob)
      return(summ)
    }
    param.B = lapply( c(unique(Subgrps),0), looper)
    param.B = do.call(rbind, param.B)
    return(param.B)
  }
  bayes_param = NULL
  for (e in unique(param.dat$estimand)){
    param.B = get_posterior_ests(param.dat=param.dat, e = e)
    bayes_param = suppressWarnings( bind_rows(bayes_param, param.B) )
  }
  bayes_param = bayes_param[order(bayes_param$Subgrps, bayes_param$estimand),]
  param.dat = suppressWarnings( 
              bind_rows( data.frame(type = "obs", param.dat),
                         data.frame(type = "bayes", bayes_param) ) )
  bayes.sim = function(mu, sd, n=100000){
    return( rnorm(n=n, mean = mu, sd = sd)  )
  }
  return( list(param.dat=param.dat, bayes.sim=bayes.sim) )
}

