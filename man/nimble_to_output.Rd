% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nimble_to_ouput.R
\name{nimble_to_output}
\alias{nimble_to_output}
\title{Run Nowcasting Model Using Nimble Package.}
\usage{
nimble_to_output(
  seed = 1,
  formula,
  data,
  model = NULL,
  family = list(),
  delay_link = NULL,
  D = NULL,
  model_window = NULL,
  aggregate = FALSE,
  nested = FALSE,
  nested.censored = NULL,
  priors = list(),
  niter,
  nburnin,
  thin,
  nchains,
  parallel = FALSE,
  formula_to_nimble
)
}
\arguments{
\item{seed}{Set seed for MCMC.}

\item{formula}{A list with formula's for model of the response variables 'delay', if you wish to fit a 'delay' model. Additionally a formula for the 'totals' response if you wish to fit a joint model.}

\item{data}{A list with data for the delay response variable, and the totals response variable (for a 'joint' model), and any covariates.}

\item{model}{Argument set to 'delay' if only modelling the delay respnse variable, or 'joint' if you wish to jointly model the delay and totals.}

\item{family}{Set distributions for model of the delay and the totals.}

\item{delay_link}{Set the link function for the delay proportions.}

\item{D}{Set the number of delays you want to model.}

\item{model_window}{Set a window length for the nowcasting model.}

\item{aggregate}{Set to TRUE if data has been aggregated (e.g. by region or age).}

\item{nested}{Set to TRUE if you with to jointly model a nested proportion of the totals.}

\item{nested.censored}{Set the number of time steps into the past the nested count of the totals is believed to be censored.}

\item{priors}{Set priors for model parameters.}

\item{niter}{Number of iterations for MCMC chains.}

\item{nburnin}{Number of burn in for MCMC chains.}

\item{thin}{Thinning for MCMC chains.}

\item{nchains}{Set number of MCMC chains.}

\item{parallel}{How you want the model to run in parallel. Default is FALSE for not running in parallel. Set to 'chains' to run over MCMC chains, or if aggregate=TRUE set to 'aggregate' to run over aggregation.}

\item{formula_to_nimble}{Function to convert model formula to nimble objects.}
}
\value{
list of nimble output and objects.
}
\description{
Run Nowcasting Model Using Nimble Package.
}
\examples{
Runs nimble code for specified nowcasting model.
}
