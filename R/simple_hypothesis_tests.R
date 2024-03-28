# This file contains tools for working with Type S/M errors under simple
# hypothosis test scenarios with normal and t distributions.
# In future, as more complex designs such as ANOVA are added, they should each
# be given their own .R file.

#' retrodesign: Calculates Power, Type S, and Type M error
#'
#' Calculates Power, Type S, and Type M error and returns them in a list or
#' df, depending on whether a single true effect size or range is provided.
#' retro_design() is faster as it uses the closed form solution from Lu et al.
#' (2018), but this function can be used for t distributions, whereas
#' retro_design() cannot. Function originally provided in Gelman and
#' Carlin (2014), modified with permission.
#'
#'
#' @param A a numeric or list, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom. df=Inf is equivalent
#' to a normal distribution.
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error.
#' @return either a list of length 3 containing the power, type s, and type M
#' error, or if A is a list, a df that is 4 by length(A), with an effect size
#' and it's corresponding power, type s, and type m errors in each row.
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(list(.2,2,20),8.1)
#' retrodesign(.5,1,df=10)
#' @export
#' @import stats
retrodesign <- function (A, s, alpha=.05, df=Inf, n.sims=10000) {
  if (s < 0){
    stop("standard errors shouldn't be negative, try again")
  }

  UseMethod("retrodesign", A)
}

#' Numeric retrodesign
#'
#' retrodesign.numeric is the S3 method of the generic retrodesign() function,
#' used when a single numeric is passed for A. Martijn Weterings kindly
#' provided code to slightly improve this in the very low N case using
#' the non-central t-distribution.
#'
#' @param A a numeric, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom. df=Inf is equivalent
#' to a normal distribution.
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @return A list of length 3 containing the power, type s, and type M
#' error.
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(2,8.1)
#' retrodesign(.5,1,df=10)
#' @export
retrodesign.numeric <- function(A, s, alpha=.05, df=Inf, n.sims=10000){
  # Just back out n/d from df; keeps all functions with the same inputs
  n = 2 / (s^2)
  d = abs(A)

  ### boundary for alpha level t-test
  tc = qt(1-alpha/2, df = df)

  ### power
  power = pt(-tc, df = df, ncp=d/sqrt(2/n)) +
    1-pt( tc, df = df, ncp=d/sqrt(2/n))

  ### s-error rate
  type_s = pt(-tc, df = df, ncp=d/sqrt(2/n))/power

  ### simulate experiments
  x0 = rnorm(n.sims,0)

  ### m-error
  type_m = sapply(d, FUN = function(di) {
    x = abs(x0+di*sqrt(n/2))
    significant = x > tc
    return(mean(x[significant == 1]/sqrt(n/2))/di)
  })
  return(list(power = power, type_s = type_s, type_m = type_m))
}

#' List retrodesign
#'
#' retrodesign.list is the S3 method of the generic retrodesign() function,
#' used when a list is passed for A.
#'
#' @param A a list, estimates of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom. df=Inf is equivalent
#' to a normal distribution.
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @return A df that is 4 by length(A), with an effect size
#' and it's corresponding power, type s, and type m errors in each row.
#' @examples
#' retrodesign(list(.2,2,20), 8.1)
#' retrodesign(list(.2,2,20), 8.1,df = 10)
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

#' retro_design: Calculates Power, Type S, and Type M error
#'
#' This function name is deprecated in favor of the more clearly named
#' retro_design_closed_form; it won't be removed in any hurry, just trying
#' to move the naming conventions to be clearer and easier to use.
#'
#' Calculates Power, Type S, and Type M error and returns them in a list or
#' df, depending on whether a single true effect size or range is provided.
#' Uses the closed form solution found for the Type-M error found by Lu et al.
#' (2018), and thus is faster than retrodesign. For t distributions, use
#' retrodesign() instead; the closed form solution only applies in the normal
#' case.
#'
#'
#' @param A a numeric or list, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @return either a list of length 3 containing the power, type s, and type M
#' error, or if A is a list, a df that is 4 by length(A), with an effect size
#' and it's corresponding power, type s, and type m errors in each row.
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(list(.2,2,20),8.1)
#' @export
#' @import stats
retro_design <- function (A, s, alpha=.05) {
  if (s < 0){
    stop("standard errors shouldn't be negative, try again")
  }

  warning("Note: this version of the function is deprecated in favor
          of the more clearly named retro_design_closed_form. We'd suggest
          switching your code to use that; note that this is the closed
          form solution")

  UseMethod("retro_design_closed_form", A)
}

#' retro_design_closed_form: Calculates Power, Type S, and Type M error
#'
#' Calculates Power, Type S, and Type M error and returns them in a list or
#' df, depending on whether a single true effect size or range is provided.
#' Uses the closed form solution found for the Type-M error found by Lu et al.
#' (2018), and thus is faster than retrodesign. For t distributions, use
#' retrodesign() instead; the closed form solution only applies in the normal
#' case.
#'
#'
#' @param A a numeric or list, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @return either a list of length 3 containing the power, type s, and type M
#' error, or if A is a list, a df that is 4 by length(A), with an effect size
#' and it's corresponding power, type s, and type m errors in each row.
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(list(.2,2,20),8.1)
#' @export
#' @import stats
retro_design_closed_form <- function (A, s, alpha=.05) {
  if (s < 0){
    stop("standard errors shouldn't be negative, try again")
  }

  UseMethod("retro_design_closed_form", A)
}

#' Numeric retro_design_closed_form
#'
#' retro_design_closed_form.numeric is the S3 method of the generic retro_design_closed_form() function,
#' used when a single numeric is passed for A.
#'
#' @param A a numeric, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @return A list of length 3 containing the power, type s, and type M
#' error.
#' @examples
#' retrodesign(1,3.28)
#' retrodesign(2,8.1)
#' @export
retro_design_closed_form.numeric <- function(A, s, alpha=.05){
  z <- qnorm(1-alpha/2)
  p.hi <- 1 - pnorm(z-A/s)
  p.lo <- pnorm(-z-A/s)
  power <- p.hi + p.lo
  typeS <- ifelse(A >= 0, p.lo/power, 1- (p.lo/power))
  lambda <- A/s

  typeM <- (dnorm(lambda + z) + dnorm(lambda - z) +
              lambda*(pnorm(lambda + z) +pnorm(lambda-z) - 1))/
    (lambda*(1 - pnorm(lambda + z) + pnorm(lambda - z)))

  return(list(power=power, type_s=typeS, type_m=abs(typeM)))
}

#' List retro_design_closed_form
#'
#' retro_design_closed_form.list is the S3 method of the generic retro_design_closed_form() function,
#' used when a list is passed for A.
#'
#' @param A a list, estimates of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @return A df that is 4 by length(A), with an effect size
#' and it's corresponding power, type s, and type m errors in each row.
#' @examples
#' retro_design(list(.2,2,20),8.1)
#' @export
retro_design_closed_form.list <- function(A, s, alpha=.05){

  list_of_lists <- lapply(A,retro_design_closed_form.numeric,s,alpha)
  matrix_form <- do.call(rbind,list_of_lists)

  mat_with_effects <- cbind(A,matrix_form)
  ret_df <- as.data.frame(mat_with_effects)

  names <- c("effect_size", "power", "type_s", "type_m")
  colnames(ret_df) <- names

  return(ret_df)
}

#' type_s
#'
#' Calculates type s error.
#'
#' @param A a numeric or list, estimate(s) of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @return either the type S, a numeric if a single A is provided, or a df
#' of length 2 by A, with the effect size and corresponding type S error in
#' each row.
#' @examples
#' type_s(1,3.28)
#' type_s(list(.2,2,20),8.1)
#' @export
#' @import stats
type_s <- function (A, s, alpha=.05) {
  if (s < 0){
    stop("standard errors shouldn't be negative, try again")
  }

  UseMethod("type_s", A)
}

#' Numeric type_s
#'
#' this is the S3 method of the generic type_s() function,
#' used when a numeric is passed for A.
#'
#' @param A a numeric, estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @return either the type S, a numeric if a single A is provided, or a df
#' of length 2 by A, with the effect size and corresponding type S error in
#' each row.
#' @examples
#' type_s(1,3.28)
#' @export
#' @import stats
type_s.numeric <- function(A, s, alpha=.05){
  z <- qt(1-alpha/2, df=Inf)
  p.hi <- 1 - pt(z-A/s, df=Inf)
  p.lo <- pt(-z-A/s, df=Inf)
  power <- p.hi + p.lo
  typeS <- ifelse(A >= 0, p.lo/power, 1- (p.lo/power))

  return(list(type_s=typeS))
}

#' List type_s
#'
#' type_s.list is the S3 method of the generic type_s() function,
#' used when a list is passed for A.
#'
#' @param A a list, estimates of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @return A df that is 2 by length(A), with an effect size
#' and it's corresponding type s errors in each row.
#' @examples
#' type_s(list(.2,2,20),8.1)
#' @export
type_s.list <- function(A, s, alpha=.05){

  list_of_lists <- lapply(A,type_s.numeric,s,alpha)
  matrix_form <- do.call(rbind,list_of_lists)

  mat_with_effects <- cbind(A,matrix_form)
  ret_df <- as.data.frame(mat_with_effects)

  names <- c("effect_size","type_s")
  colnames(ret_df) <- names

  return(ret_df)
}

#' type_m
#'
#' Calculates type m error. Is calculated using simulation, and thus supports
#' t distributions through the df parameter.
#'
#' @param A a numeric or list, estimate(s) of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the number of degrees of freedom. df=Inf is equivalent
#' to a normal distribution.
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @return either the type m error, a numeric if a single A is provided, or a df
#' of length 2 by A, with the effect size and corresponding type m error in
#' each row.
#' @examples
#' type_m(1,3.28)
#' type_m(list(.2,2,20),8.1)
#' @export
#' @import stats
type_m <- function (A, s, alpha=.05, df=Inf, n.sims=10000) {
  if (s < 0){
    stop("standard errors shouldn't be negative, try again")
  }

  UseMethod("type_m", A)
}

#' Numeric type_m
#'
#' this is the S3 method of the generic type_m() function,
#' used when a numeric is passed for A.
#'
#' @param A a numeric, estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the number of degrees of freedom. df=Inf is equivalent
#' to a normal distribution.
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @return either the type m, a numeric if a single A is provided, or a df
#' of length 2 by A, with the effect size and corresponding type m error in
#' each row.
#' @examples
#' type_m(1,3.28)
#' @export
#' @import stats
type_m.numeric <- function(A, s, alpha=.05, df=Inf, n.sims=100000){
  # Just back out n/d from df; keeps all functions with the same inputs
  n = 2 / (s^2)
  d = abs(A)

  ### boundary for alpha level t-test
  tc = qt(1-alpha/2, df = df)

  ### power
  #power = pt(-tc, df = df, ncp=d/sqrt(2/n)) +
  #  1-pt( tc, df = df, ncp=d/sqrt(2/n))

  ### s-error rate
  #type_s = pt(-tc, df = df, ncp=d/sqrt(2/n))/power

  ### simulate experiments
  x0 = rnorm(n.sims,0)

  ### m-error
  type_m = sapply(d, FUN = function(di) {
    x = abs(x0+di*sqrt(n/2))
    significant = x > tc
    return(mean(x[significant == 1]/sqrt(n/2))/di)
  })
  return(list(type_m = type_m))
}

#' List type_m
#'
#' type_m.list is the S3 method of the generic type_m() function,
#' used when a list is passed for A.
#'
#' @param A a list, estimates of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the number of degrees of freedom. df=Inf is equivalent
#' to a normal distribution.
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
#' @return A df that is 2 by length(A), with an effect size
#' and it's corresponding type m errors in each row.
#' @examples
#' type_s(list(.2,2,20),8.1)
#' @export
type_m.list <- function(A, s, alpha=.05, df=Inf, n.sims=10000){

  list_of_lists <- lapply(A,type_m.numeric,s,alpha)
  matrix_form <- do.call(rbind,list_of_lists)

  mat_with_effects <- cbind(A,matrix_form)
  ret_df <- as.data.frame(mat_with_effects)

  names <- c("effect_size","type_m")
  colnames(ret_df) <- names

  return(ret_df)
}
