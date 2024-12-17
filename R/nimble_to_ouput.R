# Define the nimble_to_output function
#' Run Nowcasting Model Using Nimble Package.
#'
#' @param seed Set seed for MCMC.
#' @param formula A list with formula's for model of the response variables 'delay', if you wish to fit a 'delay' model. Additionally a formula for the 'totals' response if you wish to fit a joint model.
#' @param data A list with data for the delay response variable, and the totals response variable (for a 'joint' model), and any covariates.
#' @param model Argument set to 'delay' if only modelling the delay respnse variable, or 'joint' if you wish to jointly model the delay and totals.
#' @param family Set distributions for model of the delay and the totals.
#' @param delay_link Set the link function for the delay proportions.
#' @param D Set the number of delays you want to model.
#' @param model_window Set a window length for the nowcasting model.
#' @param aggregate Set to TRUE if data has been aggregated (e.g. by region or age).
#' @param nested Set to TRUE if you with to jointly model a nested proportion of the totals.
#' @param nested.censored Set the number of time steps into the past the nested count of the totals is believed to be censored.
#' @param priors Set priors for model parameters.
#' @param nchains Set number of MCMC chains.
#' @param niter Number of iterations for MCMC chains.
#' @param nburnin Number of burn in for MCMC chains.
#' @param thin Thinning for MCMC chains.
#' @param parallel How you want the model to run in parallel. Default is FALSE for not running in parallel. Set to 'chains' to run over MCMC chains, or if aggregate=TRUE set to 'aggregate' to run over aggregation.
#' @param formula_to_nimble Function to convert model formula to nimble objects.
#'
#' @return list of nimble output and objects.
#' @export
#' @import nimble
#' @importFrom nimble getNimbleOption returnType
#' @importFrom stats rbeta rbinom
#'
#' @examples
#' Runs nimble code for specified nowcasting model.
nimble_to_output <- function(seed=1, formula, data, model=NULL, family=list(),
                             delay_link=NULL,  D=NULL, model_window=NULL, aggregate=FALSE,
                             nested=FALSE, nested.censored=NULL, priors=list(),
                             niter, nburnin, thin, nchains, parallel=FALSE,
                             formula_to_nimble) {
  # Load libraries.
  # library(nimble)
  #requireNamespace("nimble", quietly = TRUE) #error
  loadNamespace("nimble")
  # library(tidyverse) # error
  # library(mgcv)
 # require(nimble, quietly = TRUE)


  # Load required functions
  if(tolower(family$delay)=="gdm"|nested){
    # Define the Beta-Binomial as a distribution for NIMBLE.
    dbetabin=nimble::nimbleFunction(run=function(x=double(0),mu=double(0),phi=double(0),size=double(0),log=integer(0)){
      phi <- min(phi,1e+04) # Hard upper limit on phi for computational stability.
      returnType(double(0))
      if(x>=0&x<=size){
        return(lgamma(size+1)+lgamma(x+mu*phi)+lgamma(size-x+(1-mu)*phi)+lgamma(phi)-
                 lgamma(size+phi)-lgamma(mu*phi)-lgamma((1-mu)*phi)-lgamma(size-x+1)-lgamma(x+1))
      }else{
        return(-Inf)
      }
    })

    rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
      phi <- min(phi,1e+04) # Hard upper limit on phi for computational stability.
      pi=rbeta(1,mu*phi,(1-mu)*phi)
      returnType(double(0))
      return(rbinom(1,size,pi))
    })

    assign('dbetabin', dbetabin, envir = .GlobalEnv)
    assign('rbetabin', rbetabin, envir = .GlobalEnv)

    # Register the Beta-Binomial as a distribution.
    registerDistributions(list(dbetabin=list(
      BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))

  }


  ############ Model response variables

  # Parse the formula$totals components
  response <- as.character(formula$totals[[2]])  # Dependent variable(s)
  # Parse the formula$delay components
  response_partial <- as.character(formula$delay[[2]])  # Dependent variable(s)

  if(nested){
    response_nested <- as.character(formula$nested[[2]])  # Dependent variable(s)
  }

  ########### Set up nimble arguments

  # Generate nimble code:

  # Generate nimble code and inputs.
  if (tolower(parallel)=="aggregate"){
    # If parallel over 'aggregate' then model each aggregated variable on each cluster.
    data_parallel<-data
    data_parallel[[response]]<-data[[response]][,seed]
    data_parallel[[response_partial]]<-data[[response_partial]][,,seed]
    if(nested){
      data_parallel[[response_nested]]<-data[[response_nested]][,seed]
    }

    nimble_code_output <- formula_to_nimble(formula=formula, data=data_parallel, model=model, family=family, delay_link=delay_link,
                                            D=D, model_window=model_window, aggregate=FALSE, nested=nested, nested.censored=nested.censored,
                                            priors=priors, nchains=nchains
    )

  }else {
    nimble_code_output <- formula_to_nimble(formula=formula, data=data, model=model, family=family, delay_link=delay_link,
                                            D=D, model_window=model_window, aggregate=aggregate, nested=nested, nested.censored=nested.censored,
                                            priors=priors, nchains=nchains
    )
  }


  # The nimble code as a string
  nimble_code_string <- nimble_code_output$nimble_code

  # Convert the string to an executable NIMBLE code object
  nimble_code <-  eval(parse(text = nimble_code_string))

  # Data for nimble.
  nimble_data <- nimble_code_output$nimble_data

  # Constants for nimble.
  nimble_constants <-  nimble_code_output$nimble_constants

  # Dimensions for nimble
  nimble_dimensions <- nimble_code_output$nimble_dimensions

  # Random initial values for nimble.
  nimble_initial_values <- nimble_code_output$nimble_initial_values

  # Build the model.
  nimble_model <- nimble::nimbleModel(code=nimble_code, constants=nimble_constants,
                                      data=nimble_data, inits=nimble_initial_values)#, dimensions = nimble_dimensions)

  # Compile the model.
  nimble_compiled_model <- nimble::compileNimble(nimble_model)

  # Set which model parameters to monitor
  nimble_monitor <- nimble_code_output$nimble_monitors

  # Set up the MCMC.
  nimble_mcmc_config <- nimble::configureMCMC(nimble_model,monitors = nimble_monitor,
                                              useConjugacy = FALSE)


  nimble_mcmc<- nimble::buildMCMC(nimble_mcmc_config)

  # Compile the MCMC.
  nimble_compiled_mcmc <- nimble::compileNimble(nimble_mcmc,project=nimble_model)

  # Run MCMC.
  nimble_output_fun<- nimble::runMCMC(nimble_compiled_mcmc,niter=niter,nburnin=nburnin,thin=thin,inits=nimble_initial_values,nchains=nchains,samplesAsCodaMCMC = TRUE)

  # list(nimble_output=nimble_output_fun, nimble_code=nimble_code)
  return(list(nimble_output=nimble_output_fun, nimble_code=nimble_code, nimble_inits=nimble_initial_values, nimble_constants=nimble_constants,
              nimble_data=nimble_data,  nimble_model=nimble_model,nimble_compiled_model=nimble_compiled_model, nimble_mcmc_config=nimble_mcmc_config))
}
