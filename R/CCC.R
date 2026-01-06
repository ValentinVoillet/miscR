#' Lin's Concordance Correlation Coefficient
#'
#' @description
#' This function calculates Lin's Concordance Correlation Coefficient (CCC) to assess
#' the agreement between two variables. It provides the CCC, accuracy (bias-shift),
#' and precision (Pearson correlation), along with their asymptotic p-values and
#' confidence intervals.
#'
#' @details
#' Original implementation by **Yunda Huang (July 2006)**.
#' Updated by **Chenchen Yu (October 2017)** to correct asymptotic bounds for accuracy.
#' Reference: JASA 97, 257-70 (2002).
#'
#' @param x Numeric vector. Target values (random or fixed).
#' @param y Numeric vector. Measures of the same length as \code{x}.
#' @param targetx Character. Indication of whether \code{x} is \code{"random"} or \code{"fixed"}.
#' Default is \code{"random"}.
#' @param alpha Numeric. Two-sided type I error level for asymptotic confidence intervals.
#' Default is \code{0.05}.
#'
#' @return A \code{list} containing:
#' \item{CCC}{Concordance Correlation Coefficient}
#' \item{CCC.pvalue}{Asymptotic p-value (inverse hyperbolic tangent transformation)}
#' \item{CCC.l}{Asymptotic lower bound of CCC}
#' \item{CCC.u}{Asymptotic upper bound of CCC}
#' \item{accuracy}{Measures shift in marginal distribution (Chi_a)}
#' \item{accuracy.pvalue}{Asymptotic p-value for accuracy (logit transformation)}
#' \item{accuracy.l}{Asymptotic lower bound of accuracy}
#' \item{accuracy.u}{Asymptotic upper bound of accuracy}
#' \item{precision}{Pearson correlation coefficient (rho)}
#' \item{precision.pvalue}{Asymptotic p-value for precision}
#' \item{precision.l}{Asymptotic lower bound of precision}
#' \item{precision.u}{Asymptotic upper bound of precision}
#'
#' @export
#'
#' @importFrom stats var cor complete.cases qnorm pnorm
CCC <- function(x, y, targetx = "random", alpha = 0.05) {

  if (length(x) != length(y)) stop("The length of x and y must be equal")

  if (all(is.na(x)) | all(is.na(y))) {
    rho_c = rho_c.pvalue = rho_c_l = rho_c_u = Chi_a = Chi_a.pvalue =
      Chi_a_l = Chi_a_u = rho = rho.pvalue = rho_l = rho_u = NA
  } else {
    ok <- stats::complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]
    n <- length(x)

    # Calculate statistics
    sdx <- sqrt((stats::var(x) + 0.0001) * (n - 1) / n)
    sdy <- sqrt((stats::var(y) + 0.0001) * (n - 1) / n)
    rho <- stats::cor(x, y, use = "complete.obs", method = "pearson")

    v2 <- (mean(x) - mean(y))^2 / (sdx * sdy)
    wbar <- sdy / sdx
    Chi_a <- 2 / (wbar + 1 / wbar + v2) # accuracy
    rho_c <- Chi_a * rho # concordance correlation coefficient

    # inverse hyperbolic tangent transformation
    Z <- 0.5 * log((1 + rho_c) / (1 - rho_c))
    meanZ <- Z

    if (targetx == "random") {
      # Asymptotic inference for CCC
      temp1 <- (1 - rho^2) * (rho_c)^2 / (1 - (rho_c)^2) / rho^2
      temp2 <- 2 * v2 * (1 - rho_c) * (rho_c)^3 / (1 - (rho_c)^2)^2 / rho
      temp3 <- v2^2 * (rho_c)^4 / 2 / (1 - (rho_c)^2)^2 / rho^2
      sdZ <- sqrt((temp1 + temp2 - temp3) / (n - 2))
      Z_l <- Z - stats::qnorm(1 - alpha / 2) * sdZ
      Z_u <- Z + stats::qnorm(1 - alpha / 2) * sdZ
      rho_c_l <- (exp(2 * Z_l) - 1) / (exp(2 * Z_l) + 1)
      rho_c_u <- (exp(2 * Z_u) - 1) / (exp(2 * Z_u) + 1)
      rho_c.pvalue <- 2 * (1 - stats::pnorm(Z, mean = 0, sd = sdZ))

      # Asymptotic inference for rho
      meanZ_rho <- 0.5 * log((1 + rho) / (1 - rho))
      sdZ_rho <- sqrt(1 / (n - 3))
      Z_l <- meanZ_rho - stats::qnorm(1 - alpha / 2) * sdZ_rho
      Z_u <- meanZ_rho + stats::qnorm(1 - alpha / 2) * sdZ_rho
      rho_l <- (exp(2 * Z_l) - 1) / (exp(2 * Z_l) + 1)
      rho_u <- (exp(2 * Z_u) - 1) / (exp(2 * Z_u) + 1)
      rho.pvalue <- 2 * (1 - stats::pnorm(meanZ_rho, mean = 0, sd = sdZ_rho))

      # Asymptotic inference for accuracy
      L <- log(Chi_a / (1 - Chi_a))
      meanL <- L
      temp1 <- (Chi_a)^2 * v2 * (wbar + 1 / wbar - 2 * rho)
      temp2 <- (Chi_a)^2 * (wbar^2 + 1 / (wbar)^2 + 2 * rho^2) / 2
      temp3 <- (1 + rho^2) * (Chi_a * v2 - 1)
      sdL <- sqrt((temp1 + temp2 + temp3) / (n - 2) / (1 - Chi_a)^2)
      L_l <- L - stats::qnorm(1 - alpha / 2) * sdL
      L_u <- L + stats::qnorm(1 - alpha / 2) * sdL
      Chi_a_l <- exp(L_l) / (exp(L_l) + 1)
      Chi_a_u <- exp(L_u) / (exp(L_u) + 1)
      Chi_a.pvalue <- 2 * (1 - stats::pnorm(L, mean = 0, sd = sdL))
    }

    if (targetx == "fixed") {
      # Placeholder for future implementation
    }
  }

  return(list(
    CCC = rho_c, CCC.pvalue = rho_c.pvalue, CCC.l = rho_c_l, CCC.u = rho_c_u,
    accuracy = Chi_a, accuracy.pvalue = Chi_a.pvalue, accuracy.l = Chi_a_l, accuracy.u = Chi_a_u,
    precision = rho, precision.pvalue = rho.pvalue, precision.l = rho_l, precision.u = rho_u
  ))
}

