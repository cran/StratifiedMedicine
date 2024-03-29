### Resampling ###
plot_resample <- function(x, target=NULL) {
  plot.dat = x$resamp_dist
  if (is.null(x$Rules)) {
    plot.dat$Rules = ifelse(plot.dat$Subgrps=="ovrl", "Overall",
                            as.character(plot.dat$Subgrps))
  }
  if (!is.null(x$Rules)) {
    rules = rbind(data.frame(Subgrps="ovrl",Rules="Overall"), x$Rules)
    plot.dat = left_join(plot.dat, rules, by="Subgrps")
  }
  ## Target? ##
  if (is.null(target)) {
    t_name <- unique(plot.dat$estimand)
    t_name <- t_name[grepl("mu", t_name) | grepl("logHR", t_name) |
                       grepl("logT", t_name)]
    t_name <- t_name[1]
  }
  if (!is.null(target)) { t_name <- target }
  plot.dat = plot.dat[plot.dat$estimand==t_name,]
  if (x$param=="cox") {
    plot.dat$est = exp(plot.dat$est)
    plot.dat$estimand = "HR(A=1 vs A=0)" 
    t_name <- "HR(A=1 vs A=0)"
  }
  res = ggplot2::ggplot(plot.dat, ggplot2::aes(est, fill=Rules)) + 
    ggplot2::geom_density(alpha=0.30) +
    ggplot2::xlab( paste("Bootstrap Estimates:", t_name)  ) +
    ggplot2::facet_wrap(~Rules) +
    ggplot2::ggtitle("Bootstrap Distribution of Overall/Subgroup Estimates")+
    ggplot2::theme(plot.title=ggplot2::element_text(size=16,face="bold"),
                   axis.text.y=ggplot2::element_blank(),
                   axis.ticks.y=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_text(face="bold"),
                   axis.title=ggplot2::element_text(size=12,face="bold"))+
    ggplot2::theme_bw()
  return(res)
}