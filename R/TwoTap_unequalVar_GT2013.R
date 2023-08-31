#' Compute likelihood for the two-tap scenario with unequal spatial variance
#'
#' @param x1m measured/sensed or hypothetical position of tap 1 for which the likelihood (p given the other parameters) is computed
#' @param x2m measured/sensed or hypothetical position of tap 2
#' @param x1 (real) position of tap 1
#' @param x2 (real) position of tap 1
#' @param sigma_s1 spatial uncertainty (given as a standard deviation) for tap 1
#' @param sigma_s2 spatial uncertainty (given as a standard deviation) for tap 2
#'
#' @description
#' The function computes likelihood values for the given trajectory.
#' The code refers to the two-tap model with equal variance (spatial uncertainty) for both taps as described in Goldreich & Tong (2013), given in **Formula A17**.
#'
#'
#' @return likelihood values (see article for details or the vignette on the two-tap scenario). If x1 and x2 are vectors, a vector of likelihood values is returned.
#'
#' @export
#'
#' @examples
#' require(rabBITS)
#'
#' #XXX Example 1: compute a single point estimate XXX
#' likelihood_2Tap_UneqVar(x1m=2, x2m=4, x1=1, x2=2, sigma_s1=1, sigma_s2=2)
#'
#' #XXX Example 2: plot likelihood function for combinations of tap 1 and 2 XXX
#' library(ggplot2)
#' x1_range <- c(0, 10) #range for taps
#' x2_range <- c(0, 10)
#'
#' x1_res <- 100 #resolution for graphs
#' x2_res <- 100
#'
#' likelihoodMat <- expand.grid(x1=seq(x1_range[1], x1_range[2], length.out=x1_res), x2=seq(x2_range[1], x2_range[2], length.out=x2_res))
#' likelihoodMat$l <- likelihood_2Tap_UneqVar(x1m=2, x2m=4, x1=likelihoodMat$x1, x2=likelihoodMat$x2, sigma_s1=1, sigma_s2=2)
#' ggplot(likelihoodMat, aes(x=x1, y=x2, fill=l)) +
#'   geom_raster() +
#'   coord_fixed() +
#'   ggtitle("likelihood")
#'
#' @references Goldreich & Tong, 2013, Frontiers in Psychology

likelihood_2Tap_UneqVar <- function(x1m, x2m, x1, x2, sigma_s1, sigma_s2) {
  p_x1m_given_x1 <- 1/sqrt(2*pi)*sigma_s1 * exp(-(x1m-x1)^2/(2*sigma_s1^2))
  p_x2m_given_x2 <- 1/sqrt(2*pi)*sigma_s2 * exp(-(x2m-x2)^2/(2*sigma_s2^2))
  Likelihood <- p_x1m_given_x1*p_x2m_given_x2
  return(Likelihood)
}


#' Compute posterior probability density for the two-tap scenario with unequal spatial variance
#'
#' @param x1m measured/sensed or hypothetical position of tap 1 for which the likelihood (p given the other parameters) is computed
#' @param x2m measured/sensed or hypothetical position of tap 2
#' @param x1 (real) position of tap 1
#' @param x2 (real) position of tap 1
#' @param time_t speed prior (in units of space per time; given as a standard deviation)
#' @param sigma_s1 spatial uncertainty (given as a standard deviation) for tap 1
#' @param sigma_s2 spatial uncertainty (given as a standard deviation) for tap 2
#' @param sigma_v speed prior (in units of space per time; given as a standard deviation)
#'
#'
#' @description
#' The function computes values that are proportional to the posterior probability values for the given trajectory.
#' The code refers to the two-tap model with unequal variances (spatial uncertainties) for both taps as described in Goldreich & Tong (2013), given in **Formula A18**.
#'
#'
#' @return values that are proportional to the posterior probability for the given trajectory. If x1 and x2 are vectors, a vector of probabilities is returned.
#' Note that the returned values are not strictly speaking probabilities because they are not normalized (see article for details or the vignette on the two-tap scenario).
#'
#' @export
#'
#' @examples
#' #XXX Example 1: compute a single point estimate XXX
#' posterior_2Tap_UneqVar(x1m=2, x2m=4, x1=2, x2=4, sigma_s1=1, sigma_s2=2, sigma_v=10, time_t=0.1)
#'
#' #XXX Example 2: plot a prior distribution for combinations of tap 1 and 2 XXX
#' library(ggplot2)
#' x1_range <- c(0, 10) #range for taps
#' x2_range <- c(0, 10)
#'
#' x1_res <- 100 #resolution for graphs
#' x2_res <- 100
#'
#' posteriorMat <- expand.grid(x1=seq(x1_range[1], x1_range[2], length.out=x1_res), x2=seq(x2_range[1], x2_range[2], length.out=x2_res))
#'
#' posteriorMat$p <- posterior_2Tap_UneqVar(x1m=2, x2m=4, x1=posteriorMat$x1, x2=posteriorMat$x2, sigma_s1=1, sigma_s2=2, sigma_v=10, time_t=0.1)
#'
#' ggplot(posteriorMat, aes(x=x1, y=x2, fill=p)) +
#'   geom_raster() +
#'   coord_fixed() +
#'   ggtitle("posterior")
#'
#' @references Goldreich & Tong, 2013, Frontiers in Psychology

posterior_2Tap_UneqVar <- function(x1m, x2m, x1, x2, time_t, sigma_s1, sigma_s2, sigma_v) {
  #Posterior <- exp(-( ((x1m-x1)^2/(2*sigma_s1^2) + (x2m-x2)^2)/(2*sigma_s2^2) + (x2-x1)^2/(2*(sigma_v*time_t)^2)))
  Posterior <- exp(-( (x1m-x1)^2/(2*sigma_s1^2) + (x2m-x2)^2/(2*sigma_s2^2) + (x2-x1)^2/(2*(sigma_v*time_t)^2) ))
  return(Posterior)
}




#' Compute parameters of the posterior distribution for the two-tap scenario with unequal spatial variance
#'
#' @param x1m measured/sensed or hypothetical position of tap 1 for which the likelihood (p given the other parameters) is computed
#' @param x2m measured/sensed or hypothetical position of tap 2
#' @param time_t speed prior (in units of space per time; given as a standard deviation)
#' @param sigma_s1 spatial uncertainty (given as a standard deviation) for tap 1
#' @param sigma_s2 spatial uncertainty (given as a standard deviation) for tap 2
#' @param sigma_v speed prior (in units of space per time; given as a standard deviation)
#'
#'
#' @description
#' The function computes parameters of the 2-dimensional posterior distrbution, namely, the posterior modes for both dimensions (corresponding to taps 1 and tap 2),
#' the variance, and correlation.
#' The code refers to the two-tap model with unequal variance (spatial uncertainties) for both taps as described in Goldreich & Tong (2013), given in **Formula A19**.
#'
#' @return a list containing the posterior modes, the variances, and the correlation.
#'
#' @export
#'
#' @examples
#' #posterior parameters for the given parameters...
#' posterior_params_2Tap_UneqVar(x1m=2, x2m=4, time_t=0.1, sigma_s=1, sigma_s2=2, sigma_v=10)
#'
#' @references Goldreich & Tong, 2013, Frontiers in Psychology

posterior_params_2Tap_UneqVar <- function(x1m, x2m, time_t, sigma_s1, sigma_s2, sigma_v) {
  x1_star <- function(x1m, x2m, time_t, sigma_s1, sigma_s2, sigma_v) {
    x1_star <- x1m * (((sigma_v*time_t)^2+sigma_s2^2) / ((sigma_v*time_t)^2+sigma_s1^2+sigma_s2^2)) + x2m * (sigma_s1^2 / ((sigma_v*time_t)^2+sigma_s1^2+sigma_s2^2))
    return(x1_star)
  }

  x2_star <- function(x1m, x2m, time_t, sigma_s1, sigma_s2, sigma_v) {
    x2_star <- x1m * (sigma_s2^2 / ((sigma_v*time_t)^2+sigma_s1^2+sigma_s2^2)) + x2m * (((sigma_v*time_t)^2+sigma_s1^2) / ((sigma_v*time_t)^2+sigma_s1^2+sigma_s2^2))
    return(x2_star)
  }

  x1_sigma_square <- function(time_t, sigma_s1, sigma_s2, sigma_v) {
    sigma_square <- sigma_s1^2 * (sigma_s2^2 + (sigma_v*time_t)^2) / (sigma_s1^2 + sigma_s2^2 + (sigma_v*time_t)^2)
    return(sigma_square)
  }

  x2_sigma_square <- function(time_t, sigma_s1, sigma_s2, sigma_v) {
    sigma_square <- sigma_s2^2 * (sigma_s1^2 + (sigma_v*time_t)^2) / (sigma_s1^2 + sigma_s2^2 + (sigma_v*time_t)^2)
    return(sigma_square)
  }

  correlation <- function(time_t, sigma_s1, sigma_s2, sigma_v) {
    distr_correlation <- (sigma_s1*sigma_s2) / sqrt( (sigma_s1^2 + (sigma_v*time_t)^2) *  (sigma_s2^2 + (sigma_v*time_t)^2))
    return(distr_correlation)
  }


  return(list(x1_star=x1_star(x1m, x2m, time_t, sigma_s1, sigma_s2, sigma_v),
              x2_star=x2_star(x1m, x2m, time_t, sigma_s1, sigma_s2, sigma_v),
              x1_sigma_square=x1_sigma_square(time_t, sigma_s1, sigma_s2, sigma_v),
              x2_sigma_square=x2_sigma_square(time_t, sigma_s1, sigma_s2, sigma_v),
              correlation=correlation(time_t, sigma_s1, sigma_s2, sigma_v)))
}





