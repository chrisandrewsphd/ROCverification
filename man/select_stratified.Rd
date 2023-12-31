% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_stratified.R
\name{select_stratified}
\alias{select_stratified}
\title{Select verification sample by stratum}
\usage{
select_stratified(probs, nperstratum = 10L, nstrata = 10L, seed)
}
\arguments{
\item{probs}{Vector of probabilities of selection.}

\item{nperstratum}{How many to select from each stratum.  If this exceeds the number in a stratum, all available will be included.}

\item{nstrata}{How many strata to use.}

\item{seed}{(optional) set random seed.}
}
\value{
data.frame of selected samples with 3 variables: index of samples selected (id), probability of selection (prob), and stratum (stratum).
}
\description{
Select verification sample by stratum
}
