#' Estimation of ROC in the presence of verification bias
#'
#' @param data data.frame
#' @param dvar name of disease variable
#' @param tvar name of scoring variable
#' @param estimators which (of 5) estimators to return
#' @param dformula formula for d given t and possibly x
#' @param vformula formula for v given t and possibly x
#' @param vvar name of selection variable
#' @param xvar name of additional model covariates (deprecated)
#' @param cc thresholds at which to estimate TPR and FPR
#'
#' @return data.frame with cc, TPR and FPR for requested estimators
#' @export
#'
#' @examples # None yet
ROCverification <- function(
    data = NULL,
    dvar = NULL, 
    tvar = NULL,
    estimators = c("cc", "fi", "msi", "ipw", "dr"),
    dformula = NULL,
    vformula = NULL,
    vvar = NULL,
    xvar = NULL,
    cc = NULL) {
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
  data$vv_ <- vv <- if (is.null(vvar)) {
    as.numeric(!is.na(data[[dd]]))
  } else {
    if (!is.character(vvar)) stop("'vvar' must be character")
    if (!(vvar %in% names(data))) stop("'vvar' must be the name of a variable in 'data'")
    data[[vvar]]
  }
  if (any(is.na(vv))) stop("'data[[vvar]]' must not have missing values")

  # if cc provided, use those (sorted, unique) thresholds
  # otherwise, use (sorted, unique) values of tvar.
  # note, sort removes NAs.
  cc <- unique(sort(
    if (!is.null(cc)) cc else data[[tvar]]))
  
  # begin data.frame to return
  retval <- data.frame(cc = cc)

  # logical matrix of inequalities
  # nrows = length(cc)
  # ncols = nrow(data)
  # i,j element TRUE iff c_i<=T_j
  ineqmat <- outer(cc, data[[tvar]], FUN = function(a, b) a <= b)
  
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
  if (any(c("fi", "msi", "dr") %in% estimators)) {
    if (is.null(dformula)) stop("'dformula' needed for fi, msi, and dr estimates")
    # possible default
    # dformula <- if (is.null(xvar)) {
    #   as.formula(sprintf("%s ~ %s", dvar, tvar))
    # } else {
    #   as.formula(sprintf(
    #     "%s ~ %s + %s", dvar, tvar,
    #     paste(xvar, sep = " + ")))
    # }
    dmod <- stats::glm(dformula, family = "binomial", data = data, subset = vv_ == 1)
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
  if (any(c("ipw", "dr") %in% estimators)) {
    if (is.null(vformula)) stop("'vformula' needed for ipw and dr estimates")
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

  if ("dr" %in% estimators) {
    retval$tprdr <- apply(
      ineqmat, 1, 
      stats::weighted.mean,
      w = (vd - (vv - vpreds) * dpreds)/vpreds)
    retval$fprdr <- apply(
      ineqmat, 1,
      stats::weighted.mean,
      w = (v1md - (vv - vpreds) * (1 - dpreds))/vpreds)
    attr(retval, "aucdr") <- trapezoid(retval$fprdr, retval$tprdr)
  }  
  
  return(retval)
}
