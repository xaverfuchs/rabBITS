% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TwoTap_unequalVar_GT2013.R
\name{likelihood_2Tap_UneqVar}
\alias{likelihood_2Tap_UneqVar}
\title{Compute likelihood for the two-tap scenario with unequal spatial variance}
\usage{
likelihood_2Tap_UneqVar(x1m, x2m, x1, x2, sigma_s1, sigma_s2)
}
\arguments{
\item{x1m}{measured/sensed or hypothetical position of tap 1 for which the likelihood (p given the other parameters) is computed}

\item{x2m}{measured/sensed or hypothetical position of tap 2}

\item{x1}{(real) position of tap 1}

\item{x2}{(real) position of tap 1}

\item{sigma_s1}{spatial uncertainty (given as a standard deviation) for tap 1}

\item{sigma_s2}{spatial uncertainty (given as a standard deviation) for tap 2}
}
\value{
likelihood values (see article for details or the vignette on the two-tap scenario). If x1 and x2 are vectors, a vector of likelihood values is returned.
}
\description{
The function computes likelihood values for the given trajectory.
The code refers to the two-tap model with equal variance (spatial uncertainty) for both taps as described in Goldreich & Tong (2013), given in \strong{Formula A17}.
}
\examples{
require(rabBITS)

#XXX Example 1: compute a single point estimate XXX
likelihood_2Tap_UneqVar(x1m=2, x2m=4, x1=1, x2=2, sigma_s1=1, sigma_s2=2)

#XXX Example 2: plot likelihood function for combinations of tap 1 and 2 XXX
library(ggplot2)
x1_range <- c(0, 10) #range for taps
x2_range <- c(0, 10)

x1_res <- 100 #resolution for graphs
x2_res <- 100

likelihoodMat <- expand.grid(x1=seq(x1_range[1], x1_range[2], length.out=x1_res), x2=seq(x2_range[1], x2_range[2], length.out=x2_res))
likelihoodMat$l <- likelihood_2Tap_UneqVar(x1m=2, x2m=4, x1=likelihoodMat$x1, x2=likelihoodMat$x2, sigma_s1=1, sigma_s2=2)
ggplot(likelihoodMat, aes(x=x1, y=x2, fill=l)) +
  geom_raster() +
  coord_fixed() +
  ggtitle("likelihood")

}
\references{
Goldreich & Tong, 2013, Frontiers in Psychology
}
