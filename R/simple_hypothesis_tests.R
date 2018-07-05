# This file contains tools for working with Type S/M errors under simple
# hypothosis test scenarios.
# In future, as more complex designs such as ANOVA are added, they should each
# be given their own .R file.

#' Calculates Power, Type S, and Type M error and returns them in a list, or
#' df, depending on whether a single true effect or range is provided.
#' Function originally provided in Gelman and Carlin (2014), reused with
#' permission under the MIT license.
#'
#'
#' @param A a numeric or list, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(2,8.1)
retrodesign <- function (A, ...) {
  UseMethod("retrodesign", A)
}

#' Calculates Power, Type S, and Type M error and returns them in a list.
#' Function originally provided in Gelman and Carlin (2014), reused with
#' permission under the MIT license.
#'
#' @param A a numeric, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(2,8.1)
retrodesign.numeric <- function(A, s, alpha=.05, df=Inf, n.sims=10000){
  z <- qt(1-alpha/2, df)
  p.hi <- 1 - pt(z-A/s, df)
  p.lo <- pt(-z-A/s, df)
  power <- p.hi + p.lo
  typeS <- p.lo/power
  estimate <- A + s*rt(n.sims,df)
  significant <- abs(estimate) > s*z
  exaggeration <- mean(abs(estimate)[significant])/A
return(list(power=power, typeS=typeS, exaggeration=exaggeration))
}

#' Calculates Power, Type S, and Type M error and returns them in a list.
#' Function originally provided in Gelman and Carlin (2014), reused with
#' permission under the MIT license.
#'
#' @param A a numeric, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(2,8.1)
retrodesign.list <- function(A, s, alpha=.05, df=Inf, n.sims=10000){

  list_of_lists <- lapply(A,retrodesign.numeric,s,alpha,df,n.sims)
  matrix_form <- do.call(rbind,list_of_lists)

  mat_with_effects <- cbind(A,matrix_form)
  ret_df <- as.data.frame(mat_with_effects)

  names <- c("effect_size", "power", "type_s", "type_m")
  colnames(ret_df) <- names

  return(ret_df)
}

