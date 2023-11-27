#' Select verification sample by stratum
#'
#' @param probs Vector of probabilities of selection.
#' @param nperstratum How many to select from each stratum.  If this exceeds the number in a stratum, all available will be included.
#' @param nstrata How many strata to use.
#' @param seed (optional) set random seed.
#'
#' @return data.frame of selected samples with 3 variables: index of samples selected (id), probability of selection (prob), and stratum (stratum).
select_stratified <- function(probs, nperstratum = 10L, nstrata = 10L, seed) {
  if (!missing(seed)) set.seed(seed)

  # vector indicating in which stratum each prob is
  stratum <- cut(
    probs,
    breaks = seq(from = 0, to = 1, length.out = nstrata + 1L),
    include.lowest = TRUE)
  
  dat <- data.frame(
    id = seq_along(probs),
    prob = probs,
    stratum = stratum)
  
  list_of_dats <- split(dat, stratum)
  
  sellist <- lapply(
    list_of_dats, FUN = function(df) {
      ind <- if (nrow(df) == 0) integer(0) # none to select
      else if (nrow(df) <= nperstratum) sample(nrow(df)) # shuffle all
      else sample(nrow(df), nperstratum) # select random
      return(df[ind, ])
    })
  sel <- do.call(rbind, sellist)
  sel <- sel[sample.int(nrow(sel)),] # reorder randomly
  
  return(sel)
}

