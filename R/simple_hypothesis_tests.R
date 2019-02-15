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
#'  Carlin (2014), reused with permission.
#'
#'
#' @param A a numeric or list, an estimate of the true effect size
#' @param s a numeric, standard error of the estimate
#' @param alpha a numeric, the statistical significance threshold
#' @param df a numeric, the degrees of freedom. df=Inf is equivalent
#' to a normal distribution.
#' @param n.sims a numeric, how many times to simulate when calculating Type M
#' error
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
#' used when a single numeric is passed for A.
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
  z <- qt(1-alpha/2, df)
  p.hi <- 1 - pt(z-A/s, df)
  p.lo <- pt(-z-A/s, df)
  power <- p.hi + p.lo
  if (A >= 0) {
    typeS <- p.lo/power
  }
  else {typeS <- 1- (p.lo/power)}
  # Error suppressed below is intentional reclying when a vector is passed
  estimate <- suppressWarnings(A + s*rt(n.sims,df))
  significant <- abs(estimate) > s*z
  exaggeration <- abs(mean(abs(estimate)[significant])/A)
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

  UseMethod("retro_design", A)
}

#' Numeric retro_design
#'
#' retro_design.numeric is the S3 method of the generic retro_design() function,
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
retro_design.numeric <- function(A, s, alpha=.05){
  z <- qt(1-alpha/2, df=Inf)
  p.hi <- 1 - pt(z-A/s, df=Inf)
  p.lo <- pt(-z-A/s, df=Inf)
  power <- p.hi + p.lo
  typeS <- p.lo/power
  lambda <- A/s

  typeM <- (dt(lambda + z, df=Inf) + dt(lambda - z, df=Inf) +
              lambda*(pt(lambda + z, df=Inf) +pt(lambda-z, df=Inf) - 1))/
                  (lambda*(1 - pt(lambda + z, df=Inf) + pt(lambda - z, df=Inf)))

  return(list(power=power, typeS=typeS, typeM=abs(typeM)))
}

#' List retro_design
#'
#' retro_design.list is the S3 method of the generic retro_design() function,
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
retro_design.list <- function(A, s, alpha=.05){

  list_of_lists <- lapply(A,retro_design.numeric,s,alpha)
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
  typeS <- p.lo/power

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
#' and it's correspondingtype s errors in each row.
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
type_m.numeric <- function(A, s, alpha=.05, df=Inf, n.sims=10000){
  z <- qt(1-alpha/2, df)
  # As in retrodesign, suppressed error below is intention recyle when a vector
  # is passes
  estimate <- suppressWarnings(A + s*rt(n.sims,df))
  significant <- abs(estimate) > s*z
  exaggeration <- mean(abs(estimate)[significant])/A
  return(list(type_m=abs(exaggeration)))
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
#' and it's correspondingtype m errors in each row.
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
