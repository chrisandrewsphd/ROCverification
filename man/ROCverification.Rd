% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROCverification.R
\name{ROCverification}
\alias{ROCverification}
\title{Estimation of ROC in the presence of verification bias}
\usage{
ROCverification(
  data = NULL,
  dvar = NULL,
  tvar = NULL,
  estimators = c("cc", "fi", "msi", "ipw", "spe"),
  dformula = NULL,
  vformula = NULL,
  vvar = NULL,
  thresh = NULL
)
}
\arguments{
\item{data}{data.frame containing the variables \code{dvar}, \code{tvar}, \code{vvar}, and others in the formulae \code{dformula} and \code{vformula}.}

\item{dvar}{name of disease variable. Currently, this variable must be coded 0/1.}

\item{tvar}{name of scoring variable to be assessed for predicting \code{dvar} with ROC/AUC.}

\item{estimators}{which estimators to return. Default is all 5. Options are
\code{"cc"} for the complete case (or naive) estimator,
\code{"fi"} for the full imputation estimator,
\code{"msi"} for the mean score imputation estimator,
\code{"ipw"} for the inverse probability weighting estimator, and
\code{"spe"} for the semiparametric efficient (double robust) estimator.}

\item{dformula}{formula for d given t (and possibly other variables).
Used in \code{glm(dformula, family = "binomial")}.}

\item{vformula}{formula for v given t (and possibly other variables).
\code{glm(vformula, family = "binomial")}}

\item{vvar}{name of variable in \code{data} indicating the
verification sample (i.e., observations that were
selected to have gold standard measured).
Currently, this variable must be coded 0/1.
If no variable name is provided (\code{NULL}, default), then the
function uses \code{!is.na(data[[dvar]])} to identify verification sample.}

\item{thresh}{thresholds at which to estimate TPR and FPR}
}
\value{
data.frame with thresh, TPR and FPR for requested estimators.
Attributes are AUCs for each requested estimator
(e.g., "auccc", "aucfi", ...).
}
\description{
Estimation of ROC in the presence of verification bias
}
\examples{
dat1 <- data.frame(
  disease = rep(c(0, 1), each = 200),
  score = runif(400, min = rep(c(0, 0.4), each = 200),
   max = rep(c(0.6, 1), each = 200)),
   v = rbinom(400, size = 1, prob = rep(c(0.3, 0.7), each = 200)))
dat1$diseasev <- ifelse(dat1$v == 1, dat1$disease, NA)
out <- ROCverification(
  dat1, dvar = "diseasev", tvar = "score",
  dformula = disease ~ qlogis(score),
  estimators = c("cc", "fi"))
str(out)
attr(out, "auccc")
attr(out, "aucfi")
plot(tprcc ~ fprcc, data = out, lty = 1, type = "l")
lines(tprfi ~ fprfi, data = out, col = 2, lty = 2)
}
\references{
Alonzo and Pepe, "Assessing accuracy of a
continuous screening test in the presence of verification bias".
Appl. Statist. (2005) 54, Part 1, pp. 173–190
}
