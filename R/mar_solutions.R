fi <- function(data, dd, tt, xx = NULL, cc = NULL) {
  # doesn't handle missing data cases yet
  
  if (is.null(cc)) cc <- unique(sort(data[[tt]]))
  form <- if (is.null(xx)) {
    stats::as.formula(sprintf("%s ~ %s", dd, tt))
  } else {
    stats::as.formula(sprintf(
      "%s ~ %s + %s", dd, tt,
      paste(xx, sep = " + ")))
  }
  
  mod <- stats::glm(form, data = data, family = "binomial")
  preds <- stats::predict(mod, type = "response")
  
  ineqmat <- outer(cc, data[[tt]], FUN = function(a, b) a <= b)
  
  tprfi <- apply(ineqmat, 1, stats::weighted.mean, w = preds)
  fprfi <- apply(ineqmat, 1, stats::weighted.mean, w = 1-preds)
  
  retval <- data.frame(
    cc, tprfi, fprfi
  )
  
  attr(retval, "auc") <- trapezoid(fprfi, tprfi)
  
  return(retval)
}



msi <- function(data, dd, tt, vv, xx = NULL, cc = NULL) {
  # doesn't handle missing data cases yet
  
  if (is.null(cc)) cc <- unique(sort(data[[tt]]))
  form <- if (is.null(xx)) {
    stats::as.formula(sprintf("%s ~ %s", dd, tt))
  } else {
    stats::as.formula(sprintf(
      "%s ~ %s + %s", dd, tt,
      paste(xx, sep = " + ")))
  }
  
  mod <- stats::glm(form, data = data, family = "binomial")
  preds <- stats::predict(mod, type = "response")
  
  obspreds <- ifelse(data[[vv]], data[[dd]], preds)
  
  # tprdenom <- sum(obspreds)
  # fprdenom <- sum(1-obspreds)
  
  ineqmat <- outer(cc, data[[tt]], FUN = function(a, b) a <= b)
  
  tprmsi <- apply(ineqmat, 1, stats::weighted.mean, w = obspreds)
  fprmsi <- apply(ineqmat, 1, stats::weighted.mean, w = 1-obspreds)
  
  retval <- data.frame(
    cc, tprmsi, fprmsi
  )
  
  return(retval)
}


ipw <- function(data, dd, tt, vv, xx = NULL, cc = NULL) {
  # doesn't handle missing data cases yet
  
  if (is.null(cc)) cc <- unique(sort(data[[tt]]))
  vform <- if (is.null(xx)) {
    stats::as.formula(sprintf("%s ~ %s", vv, tt))
  } else {
    stats::as.formula(sprintf(
      "%s ~ %s + %s", vv, tt,
      paste(xx, sep = " + ")))
  }
  
  vmod <- stats::glm(vform, data = data, family = "binomial")
  vpreds <- stats::predict(vmod, type = "response")
  
  ineqmat <- outer(cc, data[[tt]], FUN = function(a, b) a <= b)
  
  tpripw <- apply(ineqmat, 1, stats::weighted.mean, w = data[[vv]]*data[[dd]]/vpreds)
  fpripw <- apply(ineqmat, 1, stats::weighted.mean, w = data[[vv]]*(1-data[[dd]])/vpreds)
  
  retval <- data.frame(
    cc, tpripw, fpripw
  )
  
  return(retval)
}


dr <- function(data, dd, tt, vv, xx = NULL, cc = NULL) {
  # doesn't handle missing data cases yet
  
  if (is.null(cc)) cc <- unique(sort(data[[tt]]))
  form <- if (is.null(xx)) {
    stats::as.formula(sprintf("%s ~ %s", dd, tt))
  } else {
    stats::as.formula(sprintf(
      "%s ~ %s + %s", dd, tt,
      paste(xx, sep = " + ")))
  }
  vform <- if (is.null(xx)) {
    stats::as.formula(sprintf("%s ~ %s", vv, tt))
  } else {
    stats::as.formula(sprintf(
      "%s ~ %s + %s", vv, tt,
      paste(xx, sep = " + ")))
  }
  
  mod <- stats::glm(form, data = data, family = "binomial")
  preds <- stats::predict(mod, type = "response")
  vmod <- stats::glm(vform, data = data, family = "binomial")
  vpreds <- stats::predict(vmod, type = "response")
  
  ineqmat <- outer(cc, data[[tt]], FUN = function(a, b) a <= b)
  
  tprdr <- apply(ineqmat, 1, stats::weighted.mean, w = data[[vv]]*data[[dd]]/vpreds - (data[[vv]]-vpreds)*preds/vpreds)
  fprdr <- apply(ineqmat, 1, stats::weighted.mean, w = data[[vv]]*(1-data[[dd]])/vpreds - (data[[vv]]-vpreds)*(1-preds)/vpreds)
  
  retval <- data.frame(
    cc, tprdr, fprdr
  )
  
  return(retval)
}

