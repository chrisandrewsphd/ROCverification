#' Compute AUC (and ROC) when verification sample selected by preliminary score,
#' Full Imputation Method (old format; consider ROCverification() instead)
#'
#' @param dat data.frame containing model variables, including scrvar
#' @param model Model to estimate selection to verification sample based on
#'   screening score
#' @param scrvar Name of screening variable within dat
#'
#' @return Area under the ROC curve estimated using Full Imputation method (see
#'   Alonzo and Pepe, 2003). Also up to 100 points of ROC curve.
#' @export
#'
#' @examples
#' dat1 <- data.frame(
#'   disease = rep(c(0, 1), each = 200),
#'   score = runif(400, min = rep(c(0, 0.6), each = 200),
#'    max = rep(c(0.4, 1), each = 200)),
#'   v = rbinom(400, size = 1, prob = rep(c(0.3, 0.7), each = 200)))
#' dat1$diseasev <- ifelse(dat1$v == 1, dat1$disease, NA)
#'
#' # rms::val.prob(dat1$score, dat1$disease)
#'
#' # rms::val.prob(dat1$score[dat1$v == 1], dat1$disease[dat1$v == 1])
#'
#' AUC_FI(dat1, model = diseasev ~ qlogis(score))
AUC_FI <- function(
    dat,
    model = disease ~ qlogis(score), 
    scrvar = "score") {
  
  mod_FI <- stats::glm(
    formula = model,
    data = dat,
    family = stats::binomial)
  
  rho_FI <- stats::predict(mod_FI, type = "response", newdata = dat)
  
  # thresholds <- sort(unique(dat[[scrvar]]))
  thresholds <- unique(stats::quantile(
    dat[[scrvar]],
    prob = seq(0, 1, length.out = 101L),
    type = 1L,
    na.rm = TRUE))

  # True Positive Rate  
  PrD1_FI <- mean(rho_FI, na.rm = TRUE)
  PrD1Tc_FI <- sapply(thresholds, FUN = function(cc) mean(rho_FI * (dat$score >= cc), na.rm = TRUE))
  TPR_FI <- PrD1Tc_FI / PrD1_FI

  # False Positive Rate
  PrD0_FI <- 1 - PrD1_FI
  PrD0Tc_FI <- sapply(thresholds, FUN = function(cc) mean((1 - rho_FI) * (dat$score >= cc), na.rm = TRUE))
  FPR_FI <- PrD0Tc_FI / PrD0_FI

  AUC_FI <- trapezoid(FPR_FI, TPR_FI)
  ROC_FI <- data.frame(threshold = thresholds, TPR_FI = TPR_FI, FPR_FI = FPR_FI)
  attr(AUC_FI, "ROC") <- ROC_FI
  
  return(AUC_FI)
}

# # model for probability of disease given screening score T
# # built using validation sample
# mod_FI <- glm(disease ~ qlogis(score), data = dat1, subset = v == 1, family = binomial)
# summary(mod_FI)
# 
# # estimate probability of disease for every person
# rho_FI <- predict(mod_FI, type = "response", newdata = dat1)
# summary(rho_FI)
# 
# # estimate prevalence: Pr(D = 1)
# PrD1_FI <- mean(rho_FI)
# 
# # joint probability of disease and large screening score, for many scores:
# #   Pr(T>=c, D = 1)
# # cs <- sort(dat1$score)
# cs <- sort(dat1$score[dat1$v == 1])
# # cs <- seq(min(dat1$score), max(dat1$score), length.out = 101L)
# PrD1Tc_FI <- sapply(cs, FUN = function(cc) mean(rho_FI * (dat1$score >= cc)))
# 
# # True Positive Rate (for many c)
# TPR_FI <- PrD1Tc_FI / PrD1_FI
# 
# plot(TPR_FI)
# 
# 
# # False Positive Rate
# PrD0_FI <- 1 - PrD1_FI
# PrD0Tc_FI <- sapply(cs, FUN = function(cc) mean((1 - rho_FI) * (dat1$score >= cc)))
# FPR_FI <- PrD0Tc_FI / PrD0_FI
# plot(FPR_FI)
# 
# plot(TPR_FI ~ FPR_FI)
# ROC_FI <- cbind(TPR_FI, FPR_FI)
# AUC_FI <- trapezoid(FPR_FI, TPR_FI)
