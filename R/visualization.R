
#' sim_plot: visualize type S/M errors
#'
#' Graphs type S/M errors resulting from a simulation using the provided
#' parameters (using the same simulation method as retrodesign()). Can
#' optionally display using ggplot.
#'
#'
#' @param A a numeric, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @param gg If TRUE and ggplot2 is installed, uses ggplot2 for graphic
#' @return A list of length 3 containing the power, type s, and type M
#' error.
#' @examples
#' sim_plot(1,3.28)
#' sim_plot(.5,1)
#' @export
#' @import stats
sim_plot <- function(A, s, alpha=.05, df=Inf, n.sims=5000,gg=TRUE){
  z <- qt(1-alpha/2, df)
  p.hi <- 1 - pt(z-A/s, df)
  p.lo <- pt(-z-A/s, df)
  power <- p.hi + p.lo
  typeS <- p.lo/power
  estimate <- A + s*rt(n.sims,df)

  if (requireNamespace("ggplot2", quietly = TRUE) & gg==TRUE){
    ggplot2::ggplot(as.data.frame(estimate),
                    ggplot2::aes(x=seq_along(estimate),y=estimate)) +
      ggplot2::geom_point(pch=ifelse((estimate>s*z) & (!estimate< -s*z),15,
                      ifelse((!estimate>s*z) & (estimate< -s*z),17,
                            ifelse((!estimate>s*z) & (!estimate< -s*z),16,16))),
              col=ifelse((!estimate>s*z) & (!estimate< -s*z),"grey","black")) +
      ggplot2::geom_hline(yintercept=s*z, size = 1) +
      ggplot2::geom_hline(yintercept=-s*z, size = 1) +
      ggplot2::geom_hline(yintercept=A,linetype = "dashed",size = 1) +
      ggplot2::xlab("n")
  }

  else {
  graphics::plot(estimate, pch=ifelse((estimate>s*z) & (!estimate< -s*z),0,
                            ifelse((!estimate>s*z) & (estimate< -s*z),2,
                            ifelse((!estimate>s*z) & (!estimate< -s*z),1,1))),
       col=ifelse((!estimate>s*z) & (!estimate< -s*z),"grey","black"))

  graphics::abline(h=s*z)
  graphics::abline(h= -s*z)
  graphics::abline(h=A,lty=3,lwd=3)
  }


  #return(estimate)
}
