# Define nimbleCast function
#' Define and Run Nowcasting Model Using Nimble
#'
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
#' @param threads If parallel=TRUE, specify number of computer cores to run over.
#' @param return.nimble_objects If TRUE all nimble objects used to run the nowcasting model will be returned in output list.
#' @param return.delay If TRUE the median posterior predictions and 95% prediction intervals for the delay response variable will be returned in output list.
#'
#' @return List of nimble output (and objects if return.nimble_objects=TRUE) and median posterior predictions of totals response variables (and delay response variable return.delay=TRUE).
#' @export
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @examples
#' Runs nowcasting model in NIMBLE and in parallel if specified.
nimbleCast <- function(formula, data, model=NULL, family=list(), delay_link=NULL,
                       D=NULL, model_window=NULL, aggregate=FALSE,
                       nested=FALSE,  nested.censored=NULL,
                       niter, nburnin, thin, nchains, priors=list(),
                       parallel=FALSE, threads=NULL,
                       return.nimble_objects=FALSE,return.delay=FALSE
) {

  require(nimble, quietly = TRUE)

  ############ Model response variables
  # Parse the formula$totalss components
  response <- as.character(formula$totals[[2]])  # Dependent variable(s)
  # Parse the formula$delay components
  response_partial <- as.character(formula$delay[[2]])  # Dependent variable(s)

  ############ Set default arguments where required
  if (is.null(model)){
    # EDIT - change default to GDM?

    if(is.null(family$totals)){
      print("As no model or family$total argument has been provided a 'delay' model is being run.")
      model <- "delay"

      if(is.null(family$delay)){
        family$delay<-"nb"
      }else if (tolower(family$delay)=="gdm"|tolower(family$delay)=="multinomial"|tolower(family$delay)=="dirichlet-multinomial"){
        print(paste0("As model argument is set to 'delay' family$delay='",family$delay,
                     "' is not compatable. Instead the argument has been set to family$delay='nb'."))
        family$delay<-"nb"
      }

      if(nested){
        stop("Argument nested=TRUE is only valid if argument model='joint'.")
      }

    }else{
      model <- "joint"
      if(is.null(family$delay)){
        print("As no family$delay has been given this has been set to 'gdm'.")
        family$delay<-"gdm"
      }
    }

    if(is.null(delay_link)){
      print("As no delay_link has been given this has been set to 'survivor.probit'.")
      delay_link<-"survivor.probit"
    }

  }else if (tolower(model)=="delay"){
    if(nested){
      stop("Argument nested=TRUE is only valid if argument model='joint'.")
    }
    if(!is.null(family$totals)){
      print("As model argument is set to 'delay' the argument family$totals is being ignored.")
      family$totals<-NULL
      # EDIT change to ignore model argument and fit joint model?
    }
    if(is.null(family$delay)){
      family$delay<-"nb"
    }else if (tolower(family$delay)=="gdm"|tolower(family$delay)=="multinomial"|tolower(family$delay)=="dirichlet-multinomial"){
      print(paste0("As model argument is set to 'delay' family$delay='",family$delay,
                   "' is not compatable. Instead the argument has been set to family$delay='nb'"))
      family$delay<-"nb"
    }
    if(is.null(delay_link)){
      print("As no delay_link has been given this has been set to 'survivor.probit'.")
      delay_link<-"survivor.probit"
    }
  }else if(tolower(model)=="joint"){
    if(is.null(family$totals)){
      print("As no family$totals has been given this has been set to 'nb'.")
      family$totals<-"nb"
    }
    if(is.null(family$delay)){
      print("As no family$delay has been given this has been set to 'gdm'.")
      family$delay<-"gdm"
    }
    if(is.null(delay_link)){
      print("As no delay_link has been given this has been set to 'survivor.probit'.")
      delay_link<-"survivor.probit"
    }
  }else(
    stop("Argument for 'model' not recognised.")
  )



  ########### Run model code

  # Check dimensions are as required for aggregate=TRUE.
  if(aggregate){
    if(!(length(dim(data[[response_partial]]))==3)){
      stop(paste("Variable aggregate set to TRUE but the partial counts response variable ", response_partial, " does not have three dimensions as required, with the third dimension representating the aggregated variable.", sep=""))
    }else{
      A<-dim(data[[response_partial]])[3]
    }
  }

  # Determine parallelisation from parallel agrument

  if(parallel==FALSE){
    # Run output function without cluster
    time_run<-system.time({
      nimble_run_output <- nimble_to_output(
        seed=1,
        formula=formula,
        data=data,
        model=model,
        family=family,
        delay_link=delay_link,
        D=D,
        model_window=model_window,
        aggregate=aggregate,
        nested=nested,
        nested.censored=nested.censored,
        priors=priors,
        niter=niter,
        nburnin=nburnin,
        thin=thin,
        nchains=nchains,
        parallel=parallel,
        formula_to_nimble=formula_to_nimble)
    })

    # Clean cluster output
    nimble_code<-nimble_run_output$nimble_code
    nimble_data<-nimble_run_output$nimble_data
    nimble_constants<-nimble_run_output$nimble_constants
    nimble_output<-nimble_run_output$nimble_output
    nimble_inits<-nimble_run_output$nimble_inits

    # Combine MCMC chains
    nimble_mcmc_output<-coda::as.mcmc.list(nimble_output)
    nimble_combined_output<-tibble::as_tibble(do.call('rbind',nimble_mcmc_output))

    # Number of posterior samples.
    n_sim <- dim(nimble_combined_output)[1]

    # Construct object to save response predictions
    if(aggregate){
      response_casts<-array(NA, dim = c(n_sim, nimble_constants$N, nimble_constants$A))
    }else{
      response_casts<-array(NA, dim = c(n_sim, nimble_constants$N))
    }

    # Nowcasts and forecasts
    if(tolower(model)=="joint"){
      # Response nowcasts
      if(aggregate){
        response_casts[,1:nimble_constants$W,] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$A))

      }else{
        response_casts[,1:nimble_constants$W] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response))%>%as.matrix()

      }
      if(nimble_constants$N>nimble_constants$W){
        # Response forecasts
        if(aggregate){
          if(tolower(family$totals)=="nb"){
            # Negative-Binomial dispersion parameters.
            theta <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('theta'))%>%as.matrix()

            # Negative-Binomial means.
            lambda <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$A))
            for(a in 1:nimble_constants$A){
              response_casts[,(nimble_constants$W+1):nimble_constants$N,a] <- stats::rnbinom(n_sim*(nimble_constants$N-nimble_constants$W), mu=lambda[,(nimble_constants$W+1):nimble_constants$N,a], size=theta[,a])
            }
          }else{
            # Poisson means.
            lambda <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$A))
            for(a in 1:nimble_constants$A){
              response_casts[,(nimble_constants$W+1):nimble_constants$N,a] <- stats::rpois(n_sim*(nimble_constants$N-nimble_constants$W), lambda[,(nimble_constants$W+1):nimble_constants$N,a])
            }
          }
        }else{
          if(tolower(family$totals)=="nb"){
            # Negative-Binomial dispersion parameters.
            theta <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('theta'))%>%as.matrix()

            # Negative-Binomial means.
            lambda <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N))
            response_casts[,(nimble_constants$W+1):nimble_constants$N] <- stats::rnbinom(n_sim*(nimble_constants$N-nimble_constants$W), mu=lambda[,(nimble_constants$W+1):nimble_constants$N], size=theta)

          }else{
            # Poisson means
            lambda <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N))
            response_casts[,(nimble_constants$W+1):nimble_constants$N] <- stats::rpois(n_sim*(nimble_constants$N-nimble_constants$W), lambda[,(nimble_constants$W+1):nimble_constants$N])

          }

        }

      }
      # Calculate response quantiles:
      if(aggregate){
        response_quantiles <- response_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','a'),value.name = response)%>%
          tidyr::spread(quantile,response)
      }else{
        response_quantiles <- response_casts%>%
          apply(c(2),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t'),value.name = response)%>%
          tidyr::spread(quantile,response)

      }

      if(return.delay==TRUE){
        if(aggregate){
          # Nowcast/forecast partial counts
          response_partial_casts <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))

          if(tolower(family$delay)=="gdm"){
            # Calculate all casts
            # GDM probabilities.
            nu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('nu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D),nimble_constants$A))
            # Dirichlet dispersion parameter
            phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D),nimble_constants$A))
            nu_phi<-array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D),nimble_constants$A))

            for(d in 1:nimble_constants$D){
              for(a in 1:nimble_constants$A){
                response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
                nu_phi[,,d,a]<-nu[,,d,a]*phi[,d,a]
                p_beta <- stats::rbeta(n_sim*sum(is.na(data[[response_partial]][,d,a])),nu_phi[,which(is.na(data[[response_partial]][,d,a])),d,a], (phi[,d,a] - nu_phi[,which(is.na(data[[response_partial]][,d,a])),d,a]))
                response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a] <- stats::rbinom(n_sim*sum(is.na(data[[response_partial]][,d,a])), p_beta, response_casts[,which(is.na(data[[response_partial]][,d,a])),d,a])
              }
            }


          }else if(tolower(family$delay)=="poisson"){
            # Calculate all casts
            # Poisson means.
            mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))
            for(d in 1:nimble_constants$D){
              for(a in 1:nimble_constants$A){
                response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
                response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu[,which(is.na(data[[response_partial]][,d,a])),d,a])
              }
            }
          }else if(tolower(family$delay)=="nb"){
            # Calculate casts
            # Negative-Binomial dispersion parameters.
            phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D,nimble_constants$A))
            # Negative-Binomial means.
            mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))
            for(d in 1:nimble_constants$D){
              for(a in 1:nimble_constants$A){
                response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
                response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu=mu[,which(is.na(data[[response_partial]][,d,a])),d,a],size=phi[,d,a])
              }
            }

          }else{
            # Sample nowcasts
            response_partial_casts <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$D+1,nimble_constants$A))

            response_partial_casts[,1:nimble_constants$W,,] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response_partial))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$D+1,nimble_constants$A))

            # calculate forcasts
            if(tolower(family$delay)=="multinomial"){
              # Multinomial probabilities.
              p <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1),nimble_constants$A))
              for(a in 1:nimble_constants$A){
                for(t in (nimble_constants$W+1):nimble_constants$N){
                  response_partial_casts[,t,1:(nimble_constants$D+1),a] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p[,t,1:(nimble_constants$D+1),a], response_casts[,t,a])
                }
              }


            }else if(tolower(family$delay)=="dirichlet-multinomial"){
              # Multinomial probabilities.
              p <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1),nimble_constants$A))
              # Dirichlet dispersion parameter
              phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D+1),nimble_constants$A))
              p_phi<-array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1),nimble_constants$A))

              for(a in 1:nimble_constants$A){
                for(t in (nimble_constants$W+1):nimble_constants$N){
                  p_phi[,t,1:(nimble_constants$D+1),a]<-p[,t,1:(nimble_constants$D+1),a]*phi[,1:(nimble_constants$D+1),a]
                  p_dir<-nimble::rdirch(n_sim*(nimble_constants$D+1),
                                        p_phi[,t,1:(nimble_constants$D+1),a])
                  response_partial_casts[,t,1:(nimble_constants$D+1),a] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p_dir, response_casts[,t,a])
                }
              }
            }

          }
        }else{
          response_partial_casts <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$D))

          # Nowcast/forecast partial counts
          if(tolower(family$delay)=="gdm"){
            # Calculate all casts
            # GDM probabilities.
            nu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('nu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D)))
            # Dirichlet dispersion parameter
            phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D)))
            nu_phi<-array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D)))

            for(d in 1:nimble_constants$D){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
              nu_phi[,,d]<-nu[,,d]*phi[,d]
              p_beta <- stats::rbeta(n_sim*sum(is.na(data[[response_partial]][,d])),nu_phi[,which(is.na(data[[response_partial]][,d])),d], (phi[,d] - nu_phi[,which(is.na(data[[response_partial]][,d])),d]))
              response_partial_casts[,which(is.na(data[[response_partial]][,d])),d] <- stats::rbinom(n_sim*sum(is.na(data[[response_partial]][,d])), p_beta, response_casts[,which(is.na(data[[response_partial]][,d])),d])
            }

          }else if(tolower(family$delay)=="poisson"){
            # Calculate casts
            # Poisson means.
            mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
            for(d in 1:nimble_constants$D){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d])),d]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d])),mu[,which(is.na(data[[response_partial]][,d])),d])

            }
          }else if(tolower(family$delay)=="nb"){
            # Calculate casts
            # Negative-Binomial dispersion parameters.
            phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D))
            # Negative-Binomial means.
            mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
            for(d in 1:nimble_constants$D){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d])),d]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d])),mu=mu[,which(is.na(data[[response_partial]][,d])),d],size=phi[,d])
            }


          }else{
            response_partial_casts <- array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))

            # Sample nowcasts
            response_partial_casts[,1:nimble_constants$W] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response_partial))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$D+1))

            # calculate forcasts
            if(tolower(family$delay)=="multinomial"){
              # Multinomial probabilities.
              p <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))
              for(t in (nimble_constants$W+1):nimble_constants$N){
                response_partial_casts[,t,1:(nimble_constants$D+1)] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p[,t,1:(nimble_constants$D+1)], response_casts[,t])
              }



            }else if(tolower(family$delay)=="dirichlet-multinomial"){
              # Multinomial probabilities.
              p <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))
              # Dirichlet dispersion parameter
              phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D+1)))
              p_phi<-array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))

              for(t in (nimble_constants$W+1):nimble_constants$N){
                p_phi[,t,1:(nimble_constants$D+1)]<-p[,t,1:(nimble_constants$D+1)]*phi[,1:(nimble_constants$D+1)]
                p_dir<-nimble::rdirch(n_sim*(nimble_constants$D+1),
                                      p_phi[,t,1:(nimble_constants$D+1)])
                response_partial_casts[,t,1:(nimble_constants$D+1)] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p_dir, response_casts[,t])
              }

            }
          }
        }
        # Calculate response partial quantiles:
        if(aggregate){
          response_partial_quantiles <- response_partial_casts%>%
            apply(c(2,3,4),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','d','a'),value.name = response_partial)%>%
            tidyr::spread(quantile,response_partial)

        }else{
          response_partial_quantiles <- response_partial_casts%>%
            apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','d','t'),value.name = response_partial)%>%
            tidyr::spread(quantile,response_partial)

        }

      }


    }else{
      # Model = delay
      # Nowcast/forecast partial counts
      # Response partial nowcasts
      if(aggregate){
        response_partial_casts <- array(NA, dim=c(n_sim, nimble_constants$N, nimble_constants$D, nimble_constants$A))
        # Calculate nowcasts/forecasts
        if(tolower(family$delay)=="nb"){
          # Negative-Binomial dispersion parameters.
          phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D,nimble_constants$A))
          # Negative-Binomial means.
          mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))
          for(d in 1:nimble_constants$D){
            for(a in 1:nimble_constants$A){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu=mu[,which(is.na(data[[response_partial]][,d,a])),d,a],size=phi[,d,a])
            }
          }

        }else{
          # Poisson means.
          mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))
          for(d in 1:nimble_constants$D){
            for(a in 1:nimble_constants$A){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu[,which(is.na(data[[response_partial]][,d,a])),d,a])
            }
          }
        }
        # Calculate response casts:
        response_casts <- apply(response_partial_casts, c(1,2,4), sum)

        # Calculate response quantiles:
        response_partial_quantiles <- response_partial_casts%>%
          apply(c(2,3,4),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','d','a'),value.name = response_partial)%>%
          tidyr::spread(quantile,response_partial)
        response_quantiles <- response_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','a'),value.name = response)%>%
          tidyr::spread(quantile,response)


      }else{
        response_partial_casts[,1:nimble_constants$W,] <- array(NA, dim=c(n_sim, nimble_constants$N, nimble_constants$D))
        # Calculate nowcasts
        if(tolower(family$delay)=="nb"){
          # Negative-Binomial dispersion parameters.
          phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D))
          # Negative-Binomial means.
          mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
          for(d in 1:nimble_constants$D){
            response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
            response_partial_casts[,which(is.na(data[[response_partial]][,d])),d]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d])),mu=mu[,which(is.na(data[[response_partial]][,d])),d],size=phi[,d])
          }

        }else{
          # Poisson means.
          mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
          for(d in 1:nimble_constants$D){
            response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
            response_partial_casts[,which(is.na(data[[response_partial]][,d])),d]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d])),mu[,which(is.na(data[[response_partial]][,d])),d])

          }
        }
        # Calculate response casts:
        response_casts <- apply(response_partial_casts, c(1,2), sum)

        # Calculate response quantiles:
        response_partial_quantiles <- response_partial_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','d'),value.name = response_partial)%>%
          tidyr::spread(quantile,response_partial)
        response_quantiles <- response_casts%>%
          apply(c(2),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t'),value.name = response)%>%
          tidyr::spread(quantile,response)
      }
    }

    # If nested calculate corrected_nested nowcasts/forecasts
    if(nested){
      # Parse the formula$nested components
      response_nested <- as.character(formula$nested[[2]])  # Dependent variable(s)
      if(aggregate){
        # Nowcast/forecast partial counts
        response_nested_casts <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$A))

        # Calculate all casts
        # GDM probabilities.
        eta<- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('eta'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$A))
        # Dirichlet dispersion parameter
        chi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('chi'))%>%as.matrix()%>%matrix(nrow=n_sim, ncol=nimble_constants$A)
        eta_chi <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$A))
        # Observed and Nowcasted nested counts
        nested_cor_index<-which(stringr::str_detect(names(nimble_combined_output),paste0(response_nested, "_corrected")))
        response_nested_casts[,1:nimble_constants$W,] <- dplyr::select(nimble_combined_output,dplyr::starts_with(paste0(response_nested, "_corrected")))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$A))
        if(nimble_constants$N>nimble_constants$W){
          for(a in 1:A){
            # Forecasted nested counts
            eta_chi[,,a]<-eta[,,a]*chi[,a]
            eta_beta <- stats::rbeta(n_sim*(nimble_constants$N-nimble_constants$W),eta_chi[,(nimble_constants$W+1):nimble_constants$N,a], (chi[,a] - eta_chi[,(nimble_constants$W+1):nimble_constants$N,a]))
            response_nested_casts[ ,(nimble_constants$W+1):nimble_constants$N,a] <- stats::rbinom(n_sim*(nimble_constants$N-nimble_constants$W), eta_beta, response_casts[,(nimble_constants$W+1):nimble_constants$N,a])
          }
        }

        # Calculate response quantiles:
        response_nested_quantiles <- response_nested_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','a'),value.name = response_partial)%>%
          tidyr::spread(quantile,response_partial)

      }else{
        response_nested_casts <- array(NA, dim=c(n_sim,nimble_constants$N))
        # Calculate all casts
        # GDM probabilities.
        eta<- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('eta'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N))
        # Dirichlet dispersion parameter
        chi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('chi'))%>%as.matrix()
        eta_chi <- array(NA, dim=c(n_sim,nimble_constants$N))
        # Observed and Nowcasted nested counts
        response_nested_casts[,1:nimble_constants$W] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(paste0(response_nested, "_corrected")))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W))
        if(nimble_constants$N>nimble_constants$W){
          #Forecasted nested counts
          eta_chi<-eta*chi
          eta_beta <- stats::rbeta(n_sim*(nimble_constants$N-nimble_constants$W),eta_chi[,(nimble_constants$W+1):nimble_constants$N], (chi - eta_chi[,(nimble_constants$W+1):nimble_constants$N]))
          response_nested_casts[ ,(nimble_constants$W+1):nimble_constants$N] <- stats::rbinom(n_sim*(nimble_constants$N-nimble_constants$W), eta_beta, response_casts[,(nimble_constants$W+1):nimble_constants$N])
        }
        response_nested_quantiles <- response_nested_casts%>%
          apply(c(2),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t'),value.name = response)%>%
          tidyr::spread(quantile,response)
      }

    }

  }else if (tolower(parallel)=="aggregate"){  # Cluster over aggregate dimension (aggregate=TRUE required).
    if(!aggregate){
      stop("Model can only be parallelised over aggregate variable when the argument 'aggregate' is TRUE. ")
    }else{
      # Cluster over number of computer cores
      # Number of available computer cores.
      n_cores<-min(parallel::detectCores(), A, threads)

      # Make Cluster for MCMC.
      cluster_package<- parallel::makeCluster(n_cores)
      # Run nimble code over cluster.
      time_run<-system.time({
        nimble_run_output <- parallel::parLapply(cl = cluster_package,
                                                 fun = nimble_to_output,
                                                 X = 1:A,
                                                 formula=formula,
                                                 data=data,
                                                 model=model,
                                                 family=family,
                                                 delay_link=delay_link,
                                                 D=D,
                                                 model_window=model_window,
                                                 aggregate=aggregate,
                                                 nested=nested,
                                                 nested.censored=nested.censored,
                                                 priors=priors,
                                                 niter=niter,
                                                 nburnin=nburnin,
                                                 thin=thin,
                                                 nchains=nchains,
                                                 parallel=parallel,
                                                 formula_to_nimble=formula_to_nimble)
      })
      parallel::stopCluster(cluster_package)

      # Clean cluster output
      nimble_code<-nimble_run_output[[1]]$nimble_code
      nimble_data<-nimble_run_output[[1]]$nimble_data
      nimble_constants<-nimble_run_output[[1]]$nimble_constants
      nimble_output<-nimble_inits<-list()
      nimble_combined_output<-nimble_mcmc_output<-list()
      for(a in 1:A){
        nimble_output[[a]]<-nimble_run_output[[a]]$nimble_output
        nimble_inits[[a]]<-nimble_run_output[[a]]$nimble_inits
        # Combine MCMC chains
        nimble_mcmc_output[[a]]<-coda::as.mcmc.list(nimble_output[[a]])
        nimble_combined_output[[a]]<-tibble::as_tibble(do.call('rbind',nimble_mcmc_output[[a]]))

      }

      # Number of posterior samples.
      n_sim <- dim(nimble_combined_output[[1]])[1]

      # Construct object to save response predictions
      response_casts<-array(NA, dim = c(n_sim, nimble_constants$N, A))

      # Nowcasts and forecasts
      if(tolower(model)=="joint"){
        # Response nowcasts
        for(a in 1:A){
          response_casts[,1:nimble_constants$W, a] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with(response))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W))
        }

        if(nimble_constants$N>nimble_constants$W){
          # Response forecasts
          if(tolower(family$totals)=="nb"){
            # Negative-Binomial dispersion parameters.
            theta<-lambda<-list()
            for(a in 1:A){
              theta[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('theta'))%>%as.matrix()
              # Negative-Binomial means.
              lambda[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N))
              response_casts[,(nimble_constants$W+1):nimble_constants$N,a] <- stats::rnbinom(n_sim*(nimble_constants$N-nimble_constants$W), mu=lambda[[a]][,(nimble_constants$W+1):nimble_constants$N], size=theta[[a]])

            }

          }else{
            # Poisson means
            lambda<-list()
            for( a in 1: A){
              lambda[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N))
              response_casts[,(nimble_constants$W+1):nimble_constants$N,a] <- stats::rpois(n_sim*(nimble_constants$N-nimble_constants$W), lambda[[a]][,(nimble_constants$W+1):nimble_constants$N])

            }


          }

        }
        response_quantiles <- response_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','a'),value.name = response)%>%
          tidyr::spread(quantile,response)

        if(return.delay==TRUE){
          # Sample nowcasts
          response_partial_casts <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$D,A))

          # Nowcast/forecast partial counts
          if(tolower(family$delay)=="gdm"){
            # Calculate all casts
            # GDM probabilities.
            nu <- phi <- nu_phi <- list()
            for(a in 1:A){
              nu[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('nu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D)))
              # Dirichlet dispersion parameter
              phi[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D)))
              nu_phi[[a]] <- array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D)))

              for(d in 1:nimble_constants$D){
                response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
                nu_phi[[a]][,,d]<-nu[[a]][,,d]*phi[[a]][,d]
                p_beta <- stats::rbeta(n_sim*sum(is.na(data[[response_partial]][,d,a])),nu_phi[[a]][,which(is.na(data[[response_partial]][,d,a])),d], (phi[[a]][,d] - nu_phi[[a]][,which(is.na(data[[response_partial]][,d,a])),d]))
                response_partial_casts[,which(is.na(data[[response_partial]][,d])),d,a] <- stats::rbinom(n_sim*sum(is.na(data[[response_partial]][,d,a])), p_beta, response_casts[,which(is.na(data[[response_partial]][,d,a])),d,a])
              }
            }

          }else if(tolower(family$delay)=="poisson"){
            # Calculate casts
            # Poisson means.
            mu <- list()
            for(a in 1:A){
              mu[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
              for(d in 1:nimble_constants$D){
                response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
                response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu[[a]][,which(is.na(data[[response_partial]][,d,a])),d])

              }
            }
          }else if(tolower(family$delay)=="nb"){
            # Calculate casts
            # Negative-Binomial dispersion parameters.
            phi <- mu <- list()
            for(a in 1:A){
              phi[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D))
              # Negative-Binomial means.
              mu[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
              for(d in 1:nimble_constants$D){
                response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
                response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu=mu[[a]][,which(is.na(data[[response_partial]][,d,a])),d],size=phi[[a]][,d])
              }
            }

          }else{

            # Sample nowcasts
            response_partial_casts <- array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1),A))

            for(a in 1:A){
              response_partial_casts[,1:nimble_constants$W,a] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with(response_partial))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$D+1))
            }

            # calculate forcasts
            if(tolower(family$delay)=="multinomial"){
              # Multinomial probabilities.
              p <- list()
              for(a in 1:A){
                p[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))
                for(t in (nimble_constants$W+1):nimble_constants$N){
                  response_partial_casts[,t,1:(nimble_constants$D+1),a] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p[[a]][,t,1:(nimble_constants$D+1)], response_casts[,t,a])
                }
              }



            }else if(tolower(family$delay)=="dirichlet-multinomial"){
              # Multinomial probabilities.
              p <- phi <- p_phi <- list()
              for(a in 1:A){
                p[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))
                # Dirichlet dispersion parameter
                phi[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D+1)))
                p_phi[[a]] <- array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))

                for(t in (nimble_constants$W+1):nimble_constants$N){
                  p_phi[[a]][,t,1:(nimble_constants$D+1)]<-p[[a]][,t,1:(nimble_constants$D+1)]*phi[[a]][,1:(nimble_constants$D+1)]
                  p_dir<-nimble::rdirch(n_sim*(nimble_constants$D+1),
                                        p_phi[[a]][,t,1:(nimble_constants$D+1)])
                  response_partial_casts[,t,1:(nimble_constants$D+1),a] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p_dir, response_casts[,t,a])
                }
              }
            }
          }
          # Calculate response partial quantiles:
          if(aggregate){
            response_partial_quantiles <- response_partial_casts%>%
              apply(c(2,3,4),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','d','a'),value.name = response_partial)%>%
              tidyr::spread(quantile,response_partial)

          }else{
            response_partial_quantiles <- response_partial_casts%>%
              apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','d','t'),value.name = response_partial)%>%
              tidyr::spread(quantile,response_partial)

          }

        }


      }else{
        # Model = delay
        # Nowcast/forecast partial counts
        # Response partial nowcasts
        response_partial_casts <- array(NA, dim=c(n_sim, nimble_constants$N, nimble_constants$D, A))
        for(a in 1:A){
          response_partial_casts[,1:nimble_constants$W,,a] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with(response))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$D))
        }
        # Calculate nowcasts
        if(tolower(family$delay)=="nb"){
          # Negative-Binomial dispersion parameters.
          phi <- mu <- list()
          for(a in 1:A){
            phi[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D))
            # Negative-Binomial means.
            mu[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
            for(d in 1:nimble_constants$D){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu=mu[[a]][,which(is.na(data[[response_partial]][,d,a])),d],size=phi[[a]][,d])
            }
          }
        }else{
          # Poisson means.
          mu <- list()
          for(a in 1:A){
            mu[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
            for(d in 1:nimble_constants$D){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d,a])), mu[[a]][,which(is.na(data[[response_partial]][,d,a])),d])

            }
          }
        }
        # Calculate response casts:
        response_casts <- apply(response_partial_casts, c(1,2,4), sum)

        # Calculate response quantiles:
        response_partial_quantiles <- response_partial_casts%>%
          apply(c(2,3,4),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','d','a'),value.name = response_partial)%>%
          tidyr::spread(quantile,response_partial)
        response_quantiles <- response_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','a'),value.name = response)%>%
          tidyr::spread(quantile,response)

      }


      if(nested){
        # Nowcast/forecast partial counts
        response_nested_casts <- array(NA, dim=c(n_sim,nimble_constants$N,A))

        # Calculate all casts
        # GDM probabilities.
        eta <- chi <- eta_chi <- list()
        for(a in 1:A){
          # Observed and Nowcasted nested counts
          response_nested_casts[,1:nimble_constants$W,a] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with(paste0(response_nested, "_corrected")))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W))
          if(nimble_constants$N>nimble_constants$W){
            # Forecast nested counts
            eta[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('eta'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N))
            # Dirichlet dispersion parameter
            chi[[a]] <- dplyr::select(as.data.frame(nimble_combined_output[[a]]),dplyr::starts_with('chi'))%>%as.matrix()
            eta_chi[[a]] <- array(NA, dim=c(n_sim,nimble_constants$N))
            eta_chi[[a]]<-eta[[a]]*chi[[a]]
            eta_beta <- stats::rbeta(n_sim*(nimble_constants$N-nimble_constants$W),eta_chi[[a]][,(nimble_constants$W+1):nimble_constants$N], (chi[[a]] - eta_chi[[a]][,(nimble_constants$W+1):nimble_constants$N]))
            response_nested_casts[ ,(nimble_constants$W+1):nimble_constants$N, a] <- stats::rbinom(n_sim*(nimble_constants$N-nimble_constants$W), eta_beta, response_casts[,(nimble_constants$W+1):nimble_constants$N])
          }
        }
        # Calculate response quantiles:
        response_nested_quantiles <- response_nested_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','a'),value.name = response_partial)%>%
          tidyr::spread(quantile,response_partial)

      }


    }
  }else if(tolower(parallel)=="chains"){  # Cluster over number of MCMC chains
    # Cluster over number of computer cores
    # Number of available computer cores.
    n_cores <- min(parallel::detectCores(), nchains, threads)

    # Make Cluster for MCMC.
    cluster_package<- parallel::makeCluster(n_cores)
    # Run nimble code over cluster.
    time_run<-system.time({
      nimble_run_output <- parallel::parLapply(cl = cluster_package,
                                               fun = nimble_to_output,
                                               X = 1:nchains,
                                               formula=formula,
                                               data=data,
                                               model=model,
                                               family=family,
                                               delay_link=delay_link,
                                               D=D,
                                               model_window=model_window,
                                               aggregate=aggregate,
                                               nested=nested,
                                               nested.censored=nested.censored,
                                               priors=priors,
                                               niter=niter,
                                               nburnin=nburnin,
                                               thin=thin,
                                               nchains=1,
                                               parallel=parallel,
                                               formula_to_nimble=formula_to_nimble)
    })
    parallel::stopCluster(cluster_package)

    # Clean cluster output
    nimble_code<-nimble_run_output[[1]]$nimble_code
    nimble_data<-nimble_run_output[[1]]$nimble_data
    nimble_constants<-nimble_run_output[[1]]$nimble_constants
    nimble_output<-list()
    nimble_inits<-list()
    for(c in 1:nchains){
      nimble_output[[c]]<-nimble_run_output[[c]]$nimble_output
      nimble_inits[[c]]<-nimble_run_output[[c]]$nimble_inits
    }
    # Combine MCMC chains
    nimble_mcmc_output<-coda::as.mcmc.list(nimble_output)
    nimble_combined_output<-tibble::as_tibble(do.call('rbind',nimble_mcmc_output))

    # Number of posterior samples.
    n_sim <- dim(nimble_combined_output)[1]

    # Construct object to save response predictions
    if(aggregate){
      response_casts<-array(NA, dim = c(n_sim, nimble_constants$N, nimble_constants$A))
    }else{
      response_casts<-matrix(NA, nrow=n_sim, ncol=nimble_constants$N)
    }

    # Nowcasts and forecasts
    if(tolower(model)=="joint"){
      # Response nowcasts
      if(aggregate){
        response_casts[,1:nimble_constants$W,] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$A))

      }else{
        response_casts[,1:nimble_constants$W] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response))%>%as.matrix()

      }
      if(nimble_constants$N>nimble_constants$W){
        # Response forecasts
        if(aggregate){
          if(tolower(family$totals)=="nb"){
            # Negative-Binomial dispersion parameters.
            theta <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('theta'))%>%as.matrix()

            # Negative-Binomial means.
            lambda <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$A))
            for(a in 1:nimble_constants$A){
              response_casts[,(nimble_constants$W+1):nimble_constants$N,a] <- stats::rnbinom(n_sim*(nimble_constants$N-nimble_constants$W), mu=lambda[,(nimble_constants$W+1):nimble_constants$N,a], size=theta[,a])
            }
          }else{
            # Poisson means.
            lambda <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$A))
            for(a in 1:nimble_constants$A){
              response_casts[,(nimble_constants$W+1):nimble_constants$N,a] <- stats::rpois(n_sim*(nimble_constants$N-nimble_constants$W), lambda[,(nimble_constants$W+1):nimble_constants$N,a])
            }
          }
        }else{
          if(tolower(family$totals)=="nb"){
            # Negative-Binomial dispersion parameters.
            theta <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('theta'))%>%as.matrix()

            # Negative-Binomial means.
            lambda <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N))
            response_casts[,(nimble_constants$W+1):nimble_constants$N] <- stats::rnbinom(n_sim*(nimble_constants$N-nimble_constants$W), mu=lambda[,(nimble_constants$W+1):nimble_constants$N], size=theta)

          }else{
            # Poisson means
            lambda <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N))
            response_casts[,(nimble_constants$W+1):nimble_constants$N] <- stats::rpois(n_sim*(nimble_constants$N-nimble_constants$W), lambda[,(nimble_constants$W+1):nimble_constants$N])

          }

        }

      }
      # Calculate response quantiles:
      if(aggregate){
        response_quantiles <- response_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','a'),value.name = response)%>%
          tidyr::spread(quantile,response)
      }else{
        response_quantiles <- response_casts%>%
          apply(c(2),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t'),value.name = response)%>%
          tidyr::spread(quantile,response)

      }

      if(return.delay==TRUE){
        if(aggregate){
          # Nowcast/forecast partial counts
          response_partial_casts <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))

          if(tolower(family$delay)=="gdm"){
            # Calculate all casts
            # GDM probabilities.
            nu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('nu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D),nimble_constants$A))
            # Dirichlet dispersion parameter
            phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D),nimble_constants$A))
            nu_phi<-array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D),nimble_constants$A))

            for(d in 1:nimble_constants$D){
              for(a in 1:nimble_constants$A){
                response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
                nu_phi[,,d,a]<-nu[,,d,a]*phi[,d,a]
                p_beta <- stats::rbeta(n_sim*sum(is.na(data[[response_partial]][,d,a])),nu_phi[,which(is.na(data[[response_partial]][,d,a])),d,a], (phi[,d,a] - nu_phi[,which(is.na(data[[response_partial]][,d,a])),d,a]))
                response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a] <- stats::rbinom(n_sim*sum(is.na(data[[response_partial]][,d,a])), p_beta, response_casts[,which(is.na(data[[response_partial]][,d,a])),d,a])
              }
            }


          }else if(tolower(family$delay)=="poisson"){
            # Calculate all casts
            # Poisson means.
            mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))
            for(d in 1:nimble_constants$D){
              for(a in 1:nimble_constants$A){
                response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
                response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu[,which(is.na(data[[response_partial]][,d,a])),d,a])
              }
            }
          }else if(tolower(family$delay)=="nb"){
            # Calculate casts
            # Negative-Binomial dispersion parameters.
            phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D,nimble_constants$A))
            # Negative-Binomial means.
            mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))
            for(d in 1:nimble_constants$D){
              for(a in 1:nimble_constants$A){
                response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
                response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu=mu[,which(is.na(data[[response_partial]][,d,a])),d,a],size=phi[,d,a])
              }
            }

          }else{
            # Sample nowcasts
            response_partial_casts <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$D+1,nimble_constants$A))

            response_partial_casts[,1:nimble_constants$W,] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response_partial))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$D+1,nimble_constants$A))

            # Calculate forecasts
            if(tolower(family$delay)=="multinomial"){
              # Multinomial probabilities.
              p <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1),nimble_constants$A))
              for(a in 1:nimble_constants$A){
                for(t in (nimble_constants$W+1):nimble_constants$N){
                  response_partial_casts[,t,1:(nimble_constants$D+1),a] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p[,t,1:(nimble_constants$D+1),a], response_casts[,t,a])
                }
              }


            }else if(tolower(family$delay)=="dirichlet-multinomial"){
              # Multinomial probabilities.
              p <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1),nimble_constants$A))
              # Dirichlet dispersion parameter
              phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D+1),nimble_constants$A))
              p_phi<-array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1),nimble_constants$A))

              for(a in 1:nimble_constants$A){
                for(t in (nimble_constants$W+1):nimble_constants$N){
                  p_phi[,t,1:(nimble_constants$D+1),a]<-p[,t,1:(nimble_constants$D+1),a]*phi[,1:(nimble_constants$D+1),a]
                  p_dir<-nimble::rdirch(n_sim*(nimble_constants$D+1),
                                        p_phi[,t,1:(nimble_constants$D+1),a])
                  response_partial_casts[,t,1:(nimble_constants$D+1),a] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p_dir, response_casts[,t,a])
                }
              }
            }

          }
        }else{
          response_partial_casts <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))

          # Nowcast/forecast partial counts
          if(tolower(family$delay)=="gdm"){
            # Calculate all casts
            # GDM probabilities.
            nu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('nu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D)))
            # Dirichlet dispersion parameter
            phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D)))
            nu_phi<-array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D)))

            for(d in 1:nimble_constants$D){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
              nu_phi[,,d]<-nu[,,d]*phi[,d]
              p_beta <- stats::rbeta(n_sim*sum(is.na(data[[response_partial]][,d])),nu_phi[,which(is.na(data[[response_partial]][,d])),d], (phi[,d] - nu_phi[,which(is.na(data[[response_partial]][,d])),d]))
              response_partial_casts[,which(is.na(data[[response_partial]][,d])),d] <- stats::rbinom(n_sim*sum(is.na(data[[response_partial]][,d])), p_beta, response_casts[,which(is.na(data[[response_partial]][,d])),d])
            }

          }else if(tolower(family$delay)=="poisson"){
            # Calculate casts
            # Poisson means.
            mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
            for(d in 1:nimble_constants$D){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d])),d]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d])),mu[,which(is.na(data[[response_partial]][,d])),d])

            }
          }else if(tolower(family$delay)=="nb"){
            # Calculate casts
            # Negative-Binomial dispersion parameters.
            phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D))
            # Negative-Binomial means.
            mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
            for(d in 1:nimble_constants$D){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d])),d]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d])),mu=mu[,which(is.na(data[[response_partial]][,d])),d],size=phi[,d])
            }


          }else{
            response_partial_casts <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))

            # Sample nowcasts
            response_partial_casts[,1:nimble_constants$W,] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response_partial))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$D+1))

            # calculate forcasts
            if(tolower(family$delay)=="multinomial"){
              # Multinomial probabilities.
              p <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))
              for(t in (nimble_constants$W+1):nimble_constants$N){
                response_partial_casts[,t,1:(nimble_constants$D+1)] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p[,t,1:(nimble_constants$D+1)], response_casts[,t])
              }



            }else if(tolower(family$delay)=="dirichlet-multinomial"){
              # Multinomial probabilities.
              p <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('p'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))
              # Dirichlet dispersion parameter
              phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,(nimble_constants$D+1)))
              p_phi<-array(NA, dim=c(n_sim,nimble_constants$N,(nimble_constants$D+1)))

              for(t in (nimble_constants$W+1):nimble_constants$N){
                p_phi[,t,1:(nimble_constants$D+1)]<-p[,t,1:(nimble_constants$D+1)]*phi[,1:(nimble_constants$D+1)]
                p_dir<-nimble::rdirch(n_sim*(nimble_constants$D+1),
                                      p_phi[,t,1:(nimble_constants$D+1)])
                response_partial_casts[,t,1:(nimble_constants$D+1)] <- nimble::rmulti(n_sim*(nimble_constants$D+1), p_dir, response_casts[,t])
              }

            }
          }
        }
        # Calculate response partial quantiles:
        if(aggregate){
          response_partial_quantiles <- response_partial_casts%>%
            apply(c(2,3,4),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','d','a'),value.name = response_partial)%>%
            tidyr::spread(quantile,response_partial)

        }else{
          response_partial_quantiles <- response_partial_casts%>%
            apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','d','t'),value.name = response_partial)%>%
            tidyr::spread(quantile,response_partial)

        }

      }


    }else{
      # Model = delay
      # Nowcast/forecast partial counts
      # Response partial nowcasts
      if(aggregate){
        response_partial_casts <- array(NA, dim=c(n_sim, nimble_constants$N, nimble_constants$D, nimble_constants$A))
        # Calculate nowcasts/forecasts
        if(tolower(family$delay)=="nb"){
          # Negative-Binomial dispersion parameters.
          phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D,nimble_constants$A))
          # Negative-Binomial means.
          mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))
          for(d in 1:nimble_constants$D){
            for(a in 1:nimble_constants$A){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu=mu[,which(is.na(data[[response_partial]][,d,a])),d,a],size=phi[,d,a])
            }
          }

        }else{
          # Poisson means.
          mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D,nimble_constants$A))
          for(d in 1:nimble_constants$D){
            for(a in 1:nimble_constants$A){
              response_partial_casts[,which(!is.na(data[[response_partial]][,d,a])),d,a]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d,a])),d,a],each=n_sim)
              response_partial_casts[,which(is.na(data[[response_partial]][,d,a])),d,a]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d,a])),mu[,which(is.na(data[[response_partial]][,d,a])),d,a])
            }
          }
        }
        # Calculate response casts:
        response_casts <- apply(response_partial_casts, c(1,2,4), sum)

        # Calculate response quantiles:
        response_partial_quantiles <- response_partial_casts%>%
          apply(c(2,3,4),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','d','a'),value.name = response_partial)%>%
          tidyr::spread(quantile,response_partial)
        response_quantiles <- response_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','a'),value.name = response)%>%
          tidyr::spread(quantile,response)


      }else{
        response_partial_casts[,nimble_constants$W,] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(response))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$D))
        # Calculate nowcasts
        if(tolower(family$delay)=="nb"){
          # Negative-Binomial dispersion parameters.
          phi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$D))
          # Negative-Binomial means.
          mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
          for(d in 1:nimble_constants$D){
            response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
            response_partial_casts[,which(is.na(data[[response_partial]][,d])),d]<- stats::rnbinom(n_sim*sum(is.na(data[[response_partial]][,d])),mu=mu[,which(is.na(data[[response_partial]][,d])),d],size=phi[,d])
          }

        }else{
          # Poisson means.
          mu <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('mu'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$D))
          for(d in 1:nimble_constants$D){
            response_partial_casts[,which(!is.na(data[[response_partial]][,d])),d]<-rep(data[[response_partial]][which(!is.na(data[[response_partial]][,d])),d],each=n_sim)
            response_partial_casts[,which(is.na(data[[response_partial]][,d])),d]<- stats::rpois(n_sim*sum(is.na(data[[response_partial]][,d])),mu[,which(is.na(data[[response_partial]][,d])),d])

          }
        }
        # Calculate response casts:
        response_casts <- apply(response_partial_casts, c(1,2), sum)

        # Calculate response quantiles:
        response_partial_quantiles <- response_partial_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','d'),value.name = response_partial)%>%
          tidyr::spread(quantile,response_partial)
        response_quantiles <- response_casts%>%
          apply(c(2),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t'),value.name = response)%>%
          tidyr::spread(quantile,response)
      }
    }

    # If nested calculate corrected_nested nowcasts/forecasts
    if(nested){
      # Parse the formula$nested components
      response_nested <- as.character(formula$nested[[2]])  # Dependent variable(s)
      if(aggregate){
        # Nowcast/forecast partial counts
        response_nested_casts <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$A))

        # Calculate all casts
        # Observed and Nowcasted nested counts
        response_nested_casts[,1:nimble_constants$W,] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(paste0(response_nested, "_corrected")))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W,nimble_constants$A))
        if(nimble_constants$N>nimble_constants$W){
          # GDM probabilities.
          eta<- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('eta'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N,nimble_constants$A))
          # Dirichlet dispersion parameter
          chi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('chi'))%>%as.matrix()%>%matrix(nrow=n_sim, ncol=nimble_constants$A)
          eta_chi <- array(NA, dim=c(n_sim,nimble_constants$N,nimble_constants$A))
          for(a in 1:A){
            # Forecasted nested counts
            eta_chi[,,a]<-eta[,,a]*chi[,a]
            eta_beta <- stats::rbeta(n_sim*(nimble_constants$N-nimble_constants$W),eta_chi[,(nimble_constants$W+1):nimble_constants$N,a], (chi[,a] - eta_chi[,(nimble_constants$W+1):nimble_constants$N,a]))
            response_nested_casts[ ,(nimble_constants$W+1):nimble_constants$N,a] <- stats::rbinom(n_sim*(nimble_constants$N-nimble_constants$W), eta_beta, response_casts[,(nimble_constants$W+1):nimble_constants$N,a])
          }
        }

        # Calculate response quantiles:
        response_nested_quantiles <- response_nested_casts%>%
          apply(c(2,3),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t','a'),value.name = response_partial)%>%
          tidyr::spread(quantile,response_partial)

      }else{
        response_nested_casts <- array(NA, dim=c(n_sim,nimble_constants$N))
        # Calculate all casts
        # Observed and Nowcasted nested counts
        response_nested_casts[,1:nimble_constants$W] <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with(paste0(response_nested, "_corrected")))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$W))
        if(nimble_constants$N>nimble_constants$W){
          # GDM probabilities.
          eta<- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('eta'))%>%as.matrix()%>%array(dim=c(n_sim,nimble_constants$N))
          # Dirichlet dispersion parameter
          chi <- dplyr::select(as.data.frame(nimble_combined_output),dplyr::starts_with('chi'))%>%as.matrix()
          eta_chi <- array(NA, dim=c(n_sim,nimble_constants$N))
          #Forecast nested counts
          eta_chi<-eta*chi
          eta_beta <- stats::rbeta(n_sim*(nimble_constants$N-nimble_constants$W),eta_chi[,(nimble_constants$W+1):nimble_constants$N], (chi - eta_chi[,(nimble_constants$W+1):nimble_constants$N]))
          response_nested_casts[ ,(nimble_constants$W+1):nimble_constants$N] <- stats::rbinom(n_sim*(nimble_constants$N-nimble_constants$W), eta_beta, response_casts[,(nimble_constants$W+1):nimble_constants$N])
        }
        response_nested_quantiles <- response_nested_casts%>%
          apply(c(2),stats::quantile,c(0.025,0.5,0.975))%>%reshape2::melt(varnames=c('quantile','t'),value.name = response)%>%
          tidyr::spread(quantile,response)
      }

    }
  }else{
    # Error message - unknown parallel argument
    stop("Argument for 'parallel' is not recognised. Please set as 'FALSE' (default), 'chains' or 'aggregate' (if aggregate=TRUE). ")


  }

  # EDIT change output to just y nowcasts/forecasts?
  # IN Functions.R -> function to extract and (if applicable) simulate posterior samples of y
  return.list<- list(output_nimble=nimble_combined_output,
                     time=paste0("Run time of nimbleCast in seconds: ", round(time_run[3])))
  # Return response casts
  response.list<-list()
  response.list[[paste0(response,"_samples")]]<-response_casts
  response.list[[paste0(response,"_95quantiles")]]<-response_quantiles
  return.list[[response]]<-response.list
  # Return partial response casts if return.delay=TRUE
  if(return.delay){
    response_partial.list<-list()
    response_partial.list[[paste0(response_partial,"_samples")]]<-response_partial_casts
    response_partial.list[[paste0(response_partial,"_95quantiles")]]<-response_partial_quantiles
    return.list[[response_partial]]<-response_partial.list
  }
  if(nested){
    response_nested.list<-list()
    response_nested.list[[paste0(response_nested, "_corrected","_samples")]]<-response_nested_casts
    response_nested.list[[paste0(response_nested, "_corrected","_95quantiles")]]<-response_nested_quantiles
    return.list[[paste0(response_nested, "_corrected")]]<-response_nested.list
  }

  # Return nimble objects if return.nimble_objects=TRUE
  if(return.nimble_objects){
    return.list$nimble_objects <- list(nimble_code=nimble_code,nimble_inits=nimble_inits,nimble_constants=nimble_constants,nimble_data=nimble_data,
                                       nimble_model=nimble_run_output[[1]]$nimble_model,nimble_compiled_model=nimble_run_output[[1]]$nimble_compiled_model, nimble_mcmc_config=nimble_run_output[[1]]$nimble_mcmc_config)
  }

  # EDIT have return.monitors options? names(nimble_output)?

  return(return.list)
}
