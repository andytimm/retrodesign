# This file contains tools for working with Type S/M errors under simple
# hypothosis test scenarios.
# In future, as more complex designs such as ANOVA are added, they should each
# be given their own .R file.

#' retrodesign: Calculates Power, Type S, and Type M error
#'
#' Calculates Power, Type S, and Type M errorand returns them in a list or
#' df, depending on whether a single true effect size or range is provided.
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
#' @return either a list of length 3 containing the power, type s, and type M
#' error, or if A is a list, a df that is 4 by length(A), with an effect size
#' and it's corresponding power, type s, and type m errors in each row.
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(list(.2,2,20),8.1)
#' @export
#' @import stats
retrodesign <- function (A, s, alpha=.05, df=Inf, n.sims=10000) {
  UseMethod("retrodesign", A)
}

#' Numeric retrodesign
#'
#' retrodesign.numeric is the S3 method of the generic retrodesign() function,
#' used when a single numeric is passed for A.
#'
#' @param A a numeric, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @return A list of length 3 containing the power, type s, and type M
#' error.
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(2,8.1)
#' @export
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

#' List retrodesign
#'
#' retrodesign.list is the S3 method of the generic retrodesign() function,
#' used when a list is passed for A.
#'
#' @param A a list, estimates of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @return A df that is 4 by length(A), with an effect size
#' and it's corresponding power, type s, and type m errors in each row.
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(2,8.1)
#' @export
retrodesign.list <- function(A, s, alpha=.05, df=Inf, n.sims=10000){

  list_of_lists <- lapply(A,retrodesign.numeric,s,alpha,df,n.sims)
  matrix_form <- do.call(rbind,list_of_lists)

  mat_with_effects <- cbind(A,matrix_form)
  ret_df <- as.data.frame(mat_with_effects)

  names <- c("effect_size", "power", "type_s", "type_m")
  colnames(ret_df) <- names

  return(ret_df)
}

