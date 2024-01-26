#' Estimation of ROC in the presence of verification bias
#'
#' @param data data.frame containing the variables \code{dvar}, \code{tvar}, \code{vvar}, and others in the formulae \code{dformula} and \code{vformula}.
#' @param dvar name of disease variable. Currently, this variable must be coded 0/1.
#' @param tvar name of scoring variable to be assessed for predicting \code{dvar} with ROC/AUC.
#' @param estimators which estimators to return. Default is all 5. Options are
#' \code{"cc"} for the complete case (or naive) estimator,
#' \code{"fi"} for the full imputation estimator,
#' \code{"msi"} for the mean score imputation estimator,
#' \code{"ipw"} for the inverse probability weighting estimator, and
#' \code{"spe"} for the semiparametric efficient (double robust) estimator.
#' @param dformula formula for d given t (and possibly other variables).
#' Used in \code{glm(dformula, family = "binomial")}.
#' @param vformula formula for v given t (and possibly other variables).
#' \code{glm(vformula, family = "binomial")} 
#' @param vvar name of variable in \code{data} indicating the 
#' verification sample (i.e., observations that were 
#' selected to have gold standard measured).
#' Currently, this variable must be coded 0/1.
#' If no variable name is provided (\code{NULL}, default), then the
#' function uses \code{!is.na(data[[dvar]])} to identify verification sample.
#' @param thresh thresholds at which to estimate TPR and FPR
#'
#' @return data.frame with thresh, TPR and FPR for requested estimators.
#' Attributes are AUCs for each requested estimator
#' (e.g., "auccc", "aucfi", ...).
#' 
#' @references Alonzo and Pepe, "Assessing accuracy of a
#'  continuous screening test in the presence of verification bias".
#'  Appl. Statist. (2005) 54, Part 1, pp. 173–190
#' @export
#'
#' @examples
#' dat1 <- data.frame(
#'   disease = rep(c(0, 1), each = 200),
#'   score = runif(400, min = rep(c(0, 0.4), each = 200),
#'    max = rep(c(0.6, 1), each = 200)),
#'    v = rbinom(400, size = 1, prob = rep(c(0.3, 0.7), each = 200)))
#' dat1$diseasev <- ifelse(dat1$v == 1, dat1$disease, NA)
#' out <- ROCverification(
#'   dat1, dvar = "diseasev", tvar = "score",
#'   dformula = disease ~ qlogis(score),
#'   estimators = c("cc", "fi"))
#' str(out)
#' attr(out, "auccc")
#' attr(out, "aucfi")
#' plot(tprcc ~ fprcc, data = out, lty = 1, type = "l")
#' lines(tprfi ~ fprfi, data = out, col = 2, lty = 2)
ROCverification <- function(
    data = NULL,
    dvar = NULL, 
    tvar = NULL,
    estimators = c("cc", "fi", "msi", "ipw", "spe"),
    dformula = NULL,
    vformula = NULL,
    vvar = NULL,
    thresh = NULL) {
  # not tested with missing data (NA) in data
  
  estimators <- match.arg(estimators, several.ok = TRUE)
  
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  
  if (is.null(tvar)) stop("'tvar' must be provided")
  if (!is.character(tvar)) stop("'tvar' must be character")
  if (!(tvar %in% names(data))) stop("'tvar' must be the name of a variable in 'data'")
  tt <- data[[tvar]]
  
  if (is.null(dvar)) stop("'dvar' must be provided")
  if (!is.character(dvar)) stop("'dvar' must be character")
  if (!(dvar %in% names(data))) stop("'dvar' must be the name of a variable in 'data'")
  dd <- data[[dvar]]
  
  # if vvar provided, use it to select those observations
  # otherwise, define vv==1 iff dd is not missing
  data$vv <- vv <- if (is.null(vvar)) { # need in both environments for scoping
    as.numeric(!is.na(dd)) # or as.numeric(!is.na(data[[dvar]]))
  } else {
    if (!is.character(vvar)) stop("'vvar' must be character")
    if (!(vvar %in% names(data))) stop("'vvar' must be the name of a variable in 'data'")
    data[[vvar]]
  }
  if (any(is.na(vv))) stop("'data[[vvar]]' must not have missing values")

  # if thresh provided, use those (sorted, unique) thresholds
  # otherwise, use (sorted, unique) values of tvar.
  # note, sort removes NAs.
  thresh <- unique(sort(
    if (!is.null(thresh)) thresh else data[[tvar]]))
  
  # begin data.frame to return
  retval <- data.frame(thresh = thresh)

  # logical matrix of inequalities
  # nrows = length(thresh)
  # ncols = nrow(data)
  # i,j element TRUE iff c_i<=T_j
  ineqmat <- outer(thresh, data[[tvar]], FUN = function(a, b) a <= b)
  
  # complete case estimator
  if ("cc" %in% estimators) {
    # handle cases where dd has NAs
    # probably should check that vv has only values 0 and 1...
    vd <- ifelse(vv == 1, dd, 0)
    v1md <- ifelse(vv == 1, 1 - dd, 0)
    retval$tprcc <- apply(ineqmat, 1, stats::weighted.mean, w = vd)
    retval$fprcc <- apply(ineqmat, 1, stats::weighted.mean, w = v1md)
    attr(retval, "auccc") <- trapezoid(retval$fprcc, retval$tprcc)
  }
  
  # estimators that require an estimate of rho_i
  if (any(c("fi", "msi", "spe") %in% estimators)) {
    if (is.null(dformula)) stop("'dformula' needed for fi, msi, and spe estimates")
    # possible default
    # dformula <- if (is.null(xvar)) {
    #   as.formula(sprintf("%s ~ %s", dvar, tvar))
    # } else {
    #   as.formula(sprintf(
    #     "%s ~ %s + %s", dvar, tvar,
    #     paste(xvar, sep = " + ")))
    # }
    dmod <- stats::glm(
      dformula, family = "binomial",
      data = data, subset = vv == 1)
    dpreds <- stats::predict(dmod, type = "response", newdata = data)
  }
  
  if ("fi" %in% estimators) {
    retval$tprfi <- apply(ineqmat, 1, stats::weighted.mean, w = dpreds)
    retval$fprfi <- apply(ineqmat, 1, stats::weighted.mean, w = 1-dpreds)
    attr(retval, "aucfi") <- trapezoid(retval$fprfi, retval$tprfi)
  }
  
  if ("msi" %in% estimators) {
    obspreds <- ifelse(vv == 1, dd, dpreds)
    retval$tprmsi <- apply(ineqmat, 1, stats::weighted.mean, w = obspreds)
    obs1mpreds <- ifelse(vv == 1, 1-dd, 1-dpreds)
    retval$fprmsi <- apply(ineqmat, 1, stats::weighted.mean, w = obs1mpreds)
    attr(retval, "aucmsi") <- trapezoid(retval$fprmsi, retval$tprmsi)
  }

  # estimators that require an estimate of pi_i
  if (any(c("ipw", "spe") %in% estimators)) {
    if (is.null(vformula)) stop("'vformula' needed for ipw and spe estimates")
    # possible default
    # vformula <- if (is.null(xvar)) {
    #   as.formula(sprintf("%s ~ %s", vvar, tvar))
    # } else {
    #   as.formula(sprintf(
    #     "%s ~ %s + %s", vvar, tvar,
    #     paste(xvar, sep = " + ")))
    # }
    vmod <- stats::glm(vformula, family = "binomial", data = data)
    vpreds <- stats::predict(vmod, type = "response")
  }
  
  if ("ipw" %in% estimators) {
    retval$tpripw <- apply(ineqmat, 1, stats::weighted.mean, w = vd/vpreds)
    retval$fpripw <- apply(ineqmat, 1, stats::weighted.mean, w = v1md/vpreds)
    attr(retval, "aucipw") <- trapezoid(retval$fpripw, retval$tpripw)
  }

  if ("spe" %in% estimators) {
    retval$tprspe <- apply(
      ineqmat, 1, 
      stats::weighted.mean,
      w = (vd - (vv - vpreds) * dpreds)/vpreds)
    retval$fprspe <- apply(
      ineqmat, 1,
      stats::weighted.mean,
      w = (v1md - (vv - vpreds) * (1 - dpreds))/vpreds)
    attr(retval, "aucspe") <- trapezoid(retval$fprspe, retval$tprspe)
  }  
  
  return(retval)
}
