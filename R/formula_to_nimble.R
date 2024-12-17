# Define the function formula_to_nimble
#' Define Nowcasting Model for Nimble.
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
#'
#' @return A list of nimble objects.
#' @export
#'
#' @examples
#' Takes model formula and gives nimble objects.
formula_to_nimble <-formula_to_nimble <- function(formula, data, model=NULL, family=list(), delay_link=list(),
                                                  D=NULL, model_window=NULL, aggregate=FALSE, nested=FALSE, nested.censored=NULL,
                                                  priors=list(), nchains=1) {

  ############  Define response variables
  # Parse the formula$totals components
  response <- as.character(formula$totals[[2]])  # Dependent variable(s)
  predictors <- as.character(formula$totals[-2])  # Independent variables

  # Clean up predictors (removes "+" and spaces)
  predictors <- stringr::str_split(predictors[2], "\\+")[[1]]
  predictors <- stringr::str_trim(predictors)

  # Parse the formula$delay components
  response_partial <- as.character(formula$delay[[2]])  # Dependent variable(s)
  predictors_partial <- as.character(formula$delay[-2])  # Independent variables

  # Clean up predictors (removes "+" and spaces)
  predictors_partial <- stringr::str_split(predictors_partial[2], "\\+")[[1]]
  predictors_partial <- stringr::str_trim(predictors_partial)

  ############ Model for totals:

  # Create lists for nimble
  nimble_constants <- list() # to store constants for nimble
  nimble_initial_values<- list() # to store initial values for nimble
  for(c in 1:nchains){
    nimble_initial_values[[c]]<-list()
  }

  #nimble_dimensions <- list() # to store dimension of predictor variables for splines
  nimble_monitors  <- list() # to store which parameters to monitor

  # # Set initial values for MCMC
  # # Set defaults:
  # inits_default<-list(
  #   # default inits for totals
  #   totals.intercept="dnorm(0, 10)",
  #   totals.dispersion="dgamma(2,0.02)",
  #   totals.spline="abs(dnorm(0,1))",
  #   totals.rw2="abs(dnorm(0,1))",
  #   totals.rw1="abs(dnorm(0,1))",
  #   totals.re="abs(dnorm(0,1))",
  #   # default inits for partial counts
  #   delay.intercept="dnorm(0, 10)",
  #   delay.dispersion="dgamma(2,0.02) ",
  #   delay.spline="abs(dnorm(0,1))",
  #   delay.rw2="abs(dnorm(0,1))",
  #   delay.rw1="abs(dnorm(0,1))",
  #   delay.re="abs(dnorm(0,1))",
  #   # default inits for nested counts
  #   nested.intercept="dnorm(0, 10)",
  #   nested.dispersion="dgamma(2,0.02)",
  #   nested.scale="dnorm(0,sd=5)",
  #   nested.spline="abs(dnorm(0,1))",
  #   nested.rw2="abs(dnorm(0,1))",
  #   nested.rw1="abs(dnorm(0,1))",
  #   nested.re="abs(dnorm(0,1))",
  #   censored.dispersion="dgamma(2,0.02)",
  #   censored.slope="dgamma(shape=10,rate=200)"
  # )
  # nimble_inits<-inits_default
  # # Replace any defaults with distributions set in inits argument:
  # nimble_inits[names(inits)]<-inits


  #use:  eval(parse(text = 9))

  # Determine if there is a forecast period in the data
  if(aggregate){
    forecast <- length(which(is.na(data[[response_partial]][,1,1])))
  }else{
    forecast <- length(which(is.na(data[[response_partial]][,1])))
  }
  # Constants for nimble
  # Model window size and temporal length:
  if(is.null(model_window)){
    nimble_constants$W <- dim(data[[response_partial]])[1] - forecast # Add constant N - maximum time steps
    nimble_constants$N <- dim(data[[response_partial]])[1] # Add constant N - maximum time steps
  }else{
    nimble_constants$W <- model_window # Add constant N - maximum time steps
    nimble_constants$N <- (model_window + forecast) # Add constant N - maximum time steps with forecast
  }

  # Delay index:
  if(is.null(D)){
    if(!(tolower(family$delay)=="gdm")){
      nimble_constants$D <- dim(data[[response_partial]])[2]-1
      D_index <- nimble_constants$D + 1

    }else{
      nimble_constants$D <- dim(data[[response_partial]])[2]  # Add constant D - maximum number of delays
      D_index <- nimble_constants$D
    }

  }else{
    if(!(tolower(family$delay)=="gdm")){
      nimble_constants$D <- (D-1)
      D_index<-nimble_constants$D+1
      Dim_D<-dim(data[[response_partial]])[2]
      if(Dim_D<D_index){
        nimble_constants$D<-Dim_D-1
        D_index<-nimble_constants$D+1
        print(paste0("Argument D is too large for dimensions of data ",response_partial, ". D has been set to ",nimble_constants$D,"."))

      }else if(Dim_D>D_index){
        if(aggregate){
          All_data<-data[[response_partial]]
          New_data<-array(NA, dim=c(nimble_constants$N, D_index, nimble_constants$A))
          New_data[,1:(D_index-1),] <- All_data[,1:(D_index-1),]
          New_data[,D_index,] <- apply(All_data[,D_index:Dim_D,],c(1,3),sum)
          data[[response_partial]] <- New_data
        }else{
          All_data<-data[[response_partial]]
          New_data<-matrix(NA, nrow = nimble_constants$N, ncol=D_index)
          New_data[,1:(D_index-1)] <- All_data[,1:(D_index-1)]
          New_data[,D_index] <- apply(All_data[,D_index:Dim_D],1,sum)
          data[[response_partial]] <- New_data
        }
      }

    }else{
      nimble_constants$D <- D  # Add constant D - maximum number of delays
      D_index<-nimble_constants$D
      Dim_D<-dim(data[[response_partial]])[2]
      if(Dim_D<D_index){
        nimble_constants$D<-Dim_D
        print(paste0("Argument D is too large for dimensions of data ",response_partial, ". D has been set to ",nimble_constants$D,"."))
      }
    }
  }

  # Start building the NIMBLE code
  nimble_code <- "nimble::nimbleCode({\n"

  # If aggregate=TRUE add additional for loop to model for 2nd dimension of response data
  # Check dimensions are as required for aggregate=TRUE.
  if(aggregate){
    if(!(length(dim(data[[response_partial]]))==3)){
      stop(paste("Variable aggregate set to TRUE but the partial counts response variable ", response_partial, " does not have three dimensions as required, with the third dimension representating the aggregated variable.", sep=""))
    }else{
      nimble_code <- paste0(nimble_code, "  for(a in 1:A) {\n")
      aggregate_space<- "    "
      aggregate_index<-", a"
      aggregate_single<-"[a]"
      nimble_constants$A<-dim(data[[response]])[2]
    }
  }else{
    if((length(dim(data[[response_partial]]))>2)){
      stop(paste("Variable aggregate set to FALSE but the partial counts response variable ", response_partial, " has more than two dimensions.", sep=""))
    }else{
      aggregate_space<- "  "
      aggregate_index<-NULL
      aggregate_single<-NULL
    }
  }

  # EDIT ADD OTHER FAMILY CHOICES ??

  # Add the likelihood for model

  if (is.null(family$total)){
    conditional=TRUE

  }else if (tolower(family$total)=="nb"){
    nimble_monitors[[response]]<-response
    nimble_monitors[["theta"]]<-"theta"
    nimble_monitors[["lambda"]]<-"lambda"

    nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:W) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space,  "  ", response, "[i", aggregate_index, "] ~ dnegbin(prob=theta", aggregate_single, "/(theta", aggregate_single, " + lambda[i", aggregate_index, "]), size=theta", aggregate_single, ")\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][["theta"]] <- abs(stats::rnorm(nimble_constants$A,0,1))

      }else{
        nimble_initial_values[[c]][["theta"]] <- abs(stats::rnorm(1,0,1))
      }
    }
    conditional=FALSE

  }else if (tolower(family$total)=="poisson"){
    nimble_monitors[[response]]<-response
    nimble_monitors[["lambda"]]<-"lambda"

    nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:W) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space,  "  ", response, "[i", aggregate_index, "] ~ dpois(lambda[i", aggregate_index, "])\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
    conditional=FALSE

  }else{
    stop("Argument family$total is not recognised.")
  }



  # Initialize the linear predictor for the totals

  nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
  linear_predictor <- paste0(aggregate_space, "  log(lambda[i", aggregate_index, "]) <- beta_0", aggregate_single, " + ")

  # List to store terms like splines, random effects, and structured terms
  spline_terms <- NULL  # to store spline term basis functions
  random_effects <- list()  # for random effects
  structured_terms <- list()  # for structured additive models
  jagam_output_totals  <- list()  # to store spline mgcv::jagam output

  # Data for nimble
  nimble_data<-list()

  jagam_output <- NULL
  for (pred in predictors) {
    if (is.na(pred)){
      print(paste0("Model with just an intercept term for the totals is being fit."))

    }else if(stringr::str_detect(pred, "s\\(")) {
      # Handle spline term with mgcv::jagam
      pred_name <- stringr::str_extract(pred, "(?<=\\().*(?=\\,)")  # extract variable inside s()

      # Use mgcv::jagam to get the spline basis for the variable
      data_spline<-data
      data_spline[response]<-NULL
      data_spline[[response]]<-stats::rnorm(nimble_constants$N,0,1)


      jagam_output <- mgcv::jagam(stats::as.formula(paste(response,"~",pred)), data = data_spline, file='blank.jags')

      # Add the spline term in the linear predictor
      spline_term_name <- paste0("f_spline_", pred_name)
      linear_predictor <- paste0(linear_predictor,"beta_spline_", pred_name, "[i", aggregate_index, "] + ", sep='')
      nimble_monitors[[paste0("beta_spline_", pred_name)]]<-paste0("beta_spline_", pred_name)

      # Add spline basis to constants.
      nimble_constants[[paste("K_beta_",pred_name,sep='')]]<-dim(jagam_output$jags.data$S1)[1]
      nimble_constants[[paste("S_beta_",pred_name,sep='')]]<-jagam_output$jags.data$S1

      # Extract the spline basis matrix.
      # Calculate using mgcv::gam incase there is missing values:
      gam_output <- mgcv::gam(stats::as.formula(paste(response,"~",pred)), data = data_spline, fit=TRUE)
      blank_gam_x<-stats::predict(gam_output, newdata=data_spline, type="lpmatrix")
      # spline_basis <- jagam_output$jags.data$X[,2:(dim(jagam_output$jags.data$S1)[1]+1)]
      spline_basis <- blank_gam_x[,2:(dim(jagam_output$jags.data$S1)[1]+1)]
      nimble_constants[[spline_term_name]]<-spline_basis

      # Add predictor to constants and define spline dimensions
      #nimble_constants[[pred_name]]<-data[[pred_name]]
      # if(aggregate){
      #   nimble_dimensions[[paste0("beta_spline_", pred_name)]]<-c(max(data[[pred_name]]),nimble_constants$A)
      #
      # }else(
      #   nimble_dimensions[[paste0("beta_spline_", pred_name)]]<-c(max(data[[pred_name]]))
      #
      # )

      # Store spline basis for NIMBLE
      spline_terms[[spline_term_name]] <- pred_name
      jagam_output_totals[[spline_term_name]] <- jagam_output

      # Add initial values for the spline coefficients
      for(c in 1:nchains){
        if(aggregate){
          nimble_initial_values[[c]][[paste0("kappa_beta_", pred_name)]] <- matrix(stats::rnorm(dim(spline_basis)[2]*nimble_constants$A, 0, 0.1), ncol=nimble_constants$A)
        }else{
          nimble_initial_values[[c]][[paste0("kappa_beta_", pred_name)]] <- stats::rnorm(dim(spline_basis)[2], 0, 0.1)
        }
      }
      nimble_monitors[[paste0("kappa_beta_", pred_name)]]<-paste0("kappa_beta_", pred_name)

    } else if (stringr::str_detect(pred, "f\\(")) {
      # Handle random effects or structured models
      pred_name <- stringr::str_extract(pred, "(?<=\\().*(?=\\,)")  # extract variable inside f()

      # Determine the type of random effect or structured model
      if (stringr::str_detect(pred, "model = \"iid\"")) {
        # Random intercept
        nimble_constants[[paste0(pred_name,"_factor")]]<-as.numeric(as.factor(data[[pred_name]]))
        nimble_constants[[paste0(pred_name,"_length")]]<-length(unique(nimble_constants[[paste0(pred_name,"_factor")]]))
        linear_predictor <- paste0(linear_predictor, "u_", pred_name, "[",pred_name,"_factor[i]", aggregate_index, "] + ")
        random_effects[[pred_name]] <- "iid"
        for(c in 1:nchains){
          if(aggregate){
            nimble_initial_values[[c]][[paste0("tau_", pred_name)]] <-  abs(stats::rnorm(nimble_constants$A,0,10))
            nimble_initial_values[[c]][[paste0("u_", pred_name)]] <- matrix((stats::rnorm(nimble_constants$A*nimble_constants$N,0,1)), ncol=nimble_constants$A)

          }else{
            nimble_initial_values[[c]][[paste0("tau_", pred_name)]] <-  abs(stats::rnorm(1,0,10))
            nimble_initial_values[[c]][[paste0("u_", pred_name)]] <- (stats::rnorm(nimble_constants$N,0,10))

          }
        }
        nimble_monitors[[paste0("tau_", pred_name)]]<-paste0("tau_", pred_name)
        nimble_monitors[[paste0("u_", pred_name)]]<-paste0("u_", pred_name)

      } else if (stringr::str_detect(pred, "model = \"rw1\"")) {
        # Random walk (RW1)
        linear_predictor <- paste0(linear_predictor, "u_rw1_", pred_name, "[i", aggregate_index, "] + ")
        structured_terms[[pred_name]] <- "rw1"
        for(c in 1:nchains){
          if(aggregate){
            nimble_initial_values[[c]][[paste0("tau_rw1_", pred_name)]] <-  abs(stats::rnorm(nimble_constants$A,0,10))
            nimble_initial_values[[c]][[paste0("u_rw1_", pred_name)]] <- matrix((stats::rnorm(nimble_constants$A*nimble_constants$N,0,1)), ncol=nimble_constants$A)

          }else{
            nimble_initial_values[[c]][[paste0("tau_rw1_", pred_name)]] <-  abs(stats::rnorm(1,0,10))
            nimble_initial_values[[c]][[paste0("u_rw1_", pred_name)]] <-  (stats::rnorm(nimble_constants$N,0,10))

          }
        }
        nimble_monitors[[paste0("tau_rw1_", pred_name)]]<-paste0("tau_rw1_", pred_name)
        nimble_monitors[[paste0("u_rw1_", pred_name)]]<-paste0("u_rw1_", pred_name)

      } else{
        stop(paste0("Argument model in ",pred," is not recognised."))
      }
    } else if (pred == "1"|pred == "I") {
      # # Intercept term
      # Included by default
      # linear_predictor <- paste0(linear_predictor, "beta_int", aggregate_single, " + ")
      # for(c in 1:nchains){
      # if(aggregate){
      #   nimble_initial_values[[c]][["beta_in"]] <- stats::rnorm(nimble_constants$A, apply(apply(nimble_data[[response_partial]],c(1,3),sum,na.rm=TRUE),2,stats::median, na.rm=TRUE) , 1)
      # }else{
      #   nimble_initial_values[[c]][["beta_int"]] <- stats::rnorm(1, stats::median(apply(nimble_data[[response_partial]], 1, sum, na.rm=TRUE), na.rm=TRUE), 1)
      # }
      # }
      # nimble_monitors[[paste0("beta_int")]]<-paste0("beta_int")

    } else {
      # Regular covariate
      linear_predictor <- paste0(linear_predictor, "beta_", pred, aggregate_single,  " * ", pred, "[i", aggregate_index, "] + ")
      nimble_data[[pred]]<-data[[pred]]
      if(is.null(data[[pred]])){
        stop(paste0("Please include covariate ", pred, " in data argument."))
      }
      for(c in 1:nchains){
        if(aggregate){
          nimble_initial_values[[c]][[paste0("beta_", pred)]] <- stats::rnorm(nimble_constants$A, 0, 0.1)
        }else{
          nimble_initial_values[[c]][[paste0("beta_", pred)]] <- stats::rnorm(1, 0, 0.1)
        }
      }
      nimble_monitors[[paste0("beta_", pred)]]<-paste0("beta_", pred)

    }
  }



  # Remove trailing " + " from the linear predictor
  linear_predictor <- substr(linear_predictor, 1, nchar(linear_predictor) - 3)
  # Add the linear predictor to the code
  nimble_code <- paste0(nimble_code, linear_predictor, "\n")
  nimble_code <- paste0(nimble_code, aggregate_space, "}\n")


  ############ Model for partial counts:

  # Add the likelihood model
  if (is.null(family$delay)){ # NOTE - family shouldn't be null if run in R_package but allows function to work by self.
    # nimble_monitors[[response_partial]] <- response_partial
    nimble_monitors[[paste0("phi")]] <- paste0("phi")
    nimble_monitors[[paste0("mu")]] <- paste0("mu")
    if(!(delay_link==tolower("log"))){
      nimble_monitors[[paste0("p")]] <- paste0("p")
    }
    family$delay <- "nb"
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][[paste0("phi")]] <- matrix(abs(stats::rnorm(D_index*nimble_constants$A, 0, 1)), ncol=nimble_constants$A)
      }else{
        nimble_initial_values[[c]][[paste0("phi")]] <- abs(stats::rnorm(D_index, 0, 1))
      }
    }

    # Loop for delay
    nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:(D+1)) {\n")
    # Loop for time
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(i in 1:N) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "  ","mu[i, d", aggregate_index, "] <- p[i, d", aggregate_index, "]*lambda[i", aggregate_index, "]\n")
    # End observed time
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
    # Loop for observed time
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(i in 1:obs_index[1", aggregate_index, "]) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "  ", response_partial, "[i, d", aggregate_index, "] ~ ")
    nimble_code <- paste0(nimble_code, "dnegbin(prob=phi[d", aggregate_index, "]/(phi[d", aggregate_index, "] + mu[i, d", aggregate_index, "]), size=phi[d", aggregate_index, "]) \n")
    # End observed time loop
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
    # End delay loop
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")


  }else if (tolower(family$delay)=="gdm"){
    # nimble_monitors[[response_partial]] <- response_partial
    nimble_monitors[[paste0("phi")]] <- paste0("phi")
    nimble_monitors[[paste0("nu")]] <- paste0("nu")

    nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:obs_index[1", aggregate_index, "]) {\n")
    # Add the likelihood for the response_partial based on the family argument
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", response_partial, "[i, 1", aggregate_index, "] ~ ")
    nimble_code <- paste0(nimble_code, "dbetabin(nu[i, 1", aggregate_index, "], phi[1", aggregate_index, "], ", response, "[i", aggregate_index, "])\n")
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][[paste0("phi")]] <- matrix(abs(stats::rnorm(D_index*nimble_constants$A, 0, 5)), ncol=nimble_constants$A)
      }else{
        nimble_initial_values[[c]][[paste0("phi")]] <- abs(stats::rnorm(D_index, 0, 5))
      }
    }
    # End time loop
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
    # Loop for delay
    nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 2:D) {\n")
    # Loop for observed time
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(i in 1:obs_index[d", aggregate_index, "]) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "  ", response_partial, "[i, d", aggregate_index, "] ~ ")
    nimble_code <- paste0(nimble_code, "dbetabin(nu[i, d", aggregate_index, "], phi[d", aggregate_index, "], ", response, "[i", aggregate_index, "]- sum(", response_partial, "[i, 1:(d-1)", aggregate_index, "]))\n")
    # End observed time loop
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
    # End delay loop
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")

  }else if(tolower(family$delay)=="nb"){
    # nimble_monitors[[response_partial]] <- response_partial
    nimble_monitors[[paste0("phi")]] <- paste0("phi")
    if(!(delay_link==tolower("log"))){
      nimble_monitors[[paste0("p")]] <- paste0("p")
    }
    nimble_monitors[[paste0("mu")]] <- paste0("mu")
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][[paste0("phi")]] <- matrix(abs(stats::rnorm(D_index*nimble_constants$A, 0, 5)), ncol=nimble_constants$A)
      }else{
        nimble_initial_values[[c]][[paste0("phi")]] <- abs(stats::rnorm(D_index, 0, 5))
      }
    }

    # Loop for delay
    nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:(D+1)) {\n")

    if (!(tolower(delay_link)=="log")){
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(i in 1:N) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "  ","mu[i, d", aggregate_index, "] <- p[i, d", aggregate_index, "]*lambda[i", aggregate_index, "]\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")

    }
    # Loop for observed time
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(i in 1:obs_index[1", aggregate_index, "]) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "  ", response_partial, "[i, d", aggregate_index, "] ~ ")
    nimble_code <- paste0(nimble_code, "dnegbin(prob=phi[d", aggregate_index, "]/(phi[d", aggregate_index, "] + mu[i, d", aggregate_index, "]), size=phi[d", aggregate_index, "]) \n")
    # End observed time loop
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
    # End delay loop
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")

  }else if(tolower(family$delay)=="multinomial"){
    nimble_monitors[[response_partial]] <- response_partial
    nimble_monitors[[paste0("p")]] <- paste0("p")

    # Loop for time
    nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:W) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", response_partial, "[i, 1:(D+1)", aggregate_index, "] ~ ")
    nimble_code <- paste0(nimble_code, "dmulti(prob=p[i, 1:(D+1)", aggregate_index, "], size=", response, "[i", aggregate_index, "]) \n")
    # End time loop
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")



  }else if(tolower(family$delay)=="dirichlet-multinomial"){
    nimble_monitors[[response_partial]] <- response_partial
    nimble_monitors[[paste0("p_dir")]] <- paste0("p_dir")
    nimble_monitors[[paste0("phi")]] <- paste0("phi")
    nimble_monitors[[paste0("p")]] <- paste0("p")
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][[paste0("phi")]] <- matrix(abs(stats::rnorm(D_index*nimble_constants$A, 0, 5)), ncol=nimble_constants$A)
        nimble_initial_values[[c]][[paste0("p_dir")]] <- array(1/D_index, dim=c(nimble_constants$N,D_index,nimble_constants$A))

      }else{
        nimble_initial_values[[c]][[paste0("phi")]] <- abs(stats::rnorm(D_index, 0, 5))
        nimble_initial_values[[c]][[paste0("p_dir")]] <- array(1/D_index, dim=c(nimble_constants$N,D_index)) # EDIT - make p initial value random? (needs to sum to 1)

      }
    }
    # Loop for time
    nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:W) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "p_phi[i, 1:(D+1)", aggregate_index, "] <- p[i, 1:(D+1)", aggregate_index,"]*phi[1:(D+1)", aggregate_index,"]\n" )
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "p_dir[i, 1:(D+1)", aggregate_index, "] ~ ddirch(p_phi[i, 1:(D+1)", aggregate_index, "])\n")
    nimble_code <- paste0(nimble_code, aggregate_space, aggregate_space, "  ", response_partial, "[i, 1:(D+1)", aggregate_index, "] ~ ")
    nimble_code <- paste0(nimble_code, "dmulti(p_dir[i, 1:(D+1)", aggregate_index, "], ", response, "[i", aggregate_index, "])\n")
    # End time loop
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")


  }else if(tolower(family$delay)=="poisson"){
    #nimble_monitors[[response_partial]] <- response_partial
    nimble_monitors[[paste0("mu")]] <- paste0("mu")
    if(!(delay_link==tolower("log"))){
      nimble_monitors[[paste0("p")]] <- paste0("p")
    }
    # Loop for delay
    nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:(D+1)) {\n")
    if (!(tolower(delay_link)=="log")){
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(i in 1:N) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "  ","mu[i, d", aggregate_index, "] <- p[i, d", aggregate_index, "]*lambda[i", aggregate_index, "]\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")

    }
    # Loop for observed time
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(i in 1:obs_index[1", aggregate_index, "]) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "  ", response_partial, "[i, d", aggregate_index, "] ~ ")
    nimble_code <- paste0(nimble_code, "dpois(mu[i, d", aggregate_index, "]) \n")
    # End observed time loop
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
    # End delay loop
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")

  }else{
    stop("Argument for family$delay not recognised.")
  }


  link_function<-NULL
  # Create linear predictor for response for the partial counts
  if(is.null(delay_link)){
    link_function<-"probit"
    # EDIT? default delay link is survivor
    if(!(tolower(family$delay)=="gdm")){
      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
      # Calculate relative proportions nu from cumulative proportions S
      nimble_code <- paste0(nimble_code, aggregate_space, "  ","p[i, 1", aggregate_index, "] <- S[i, 1", aggregate_index, "]\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 2:D) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "    ","p[i, d", aggregate_index, "] <- (S[i, d", aggregate_index, "] - S[i, d-1", aggregate_index, "])\n")
      # End delay loop
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ","p[i, D+1", aggregate_index, "] <- 1 - S[i, D", aggregate_index, "]\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:D) {\n")
      # Intialise linear predictor for delay
      linear_predictor_delay<-paste0(aggregate_space, "  ", "    ","probit(S[i, d", aggregate_index, "]) <- alpha_0[d", aggregate_index, "] + ")

    }else{
      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
      # Calculate relative proportions nu from cumulative proportions S
      nimble_code <- paste0(nimble_code, aggregate_space, "  ","nu[i, 1", aggregate_index, "] <- S[i, 1", aggregate_index, "]\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 2:D) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "    ","nu[i, d", aggregate_index, "] <- (S[i, d", aggregate_index, "] - S[i, d-1", aggregate_index, "])/(1 - S[i, d-1", aggregate_index, "])\n")
      # End delay loop
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:D) {\n")
      # Intialise linear predictor for delay
      linear_predictor_delay<-paste0(aggregate_space, "  ", "    ","probit(S[i, d", aggregate_index, "]) <- alpha_0[d", aggregate_index, "] + ")

    }

  } else if(stringr::str_detect(tolower(delay_link), "survivor")) {
    if(nimble_constants$D<2){
      stop(paste0("Please increase the argument for D such that it is greater than 2 or choose an alternative link function to 'survivor'."))
    }
    link_function <- stringr::str_extract(delay_link, '\\b\\w+$') # extract variable after "-"
    if(link_function=="survivor"){
      link_function<-"probit"
    }
    if(!(tolower(family$delay)=="gdm")){
      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
      # Calculate relative proportions nu from cumulative proportions S
      nimble_code <- paste0(nimble_code, aggregate_space, "  ","p[i, 1", aggregate_index, "] <- S[i, 1", aggregate_index, "]\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 2:D) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "    ","p[i, d", aggregate_index, "] <- (S[i, d", aggregate_index, "] - S[i, d-1", aggregate_index, "])\n")
      # End delay loop
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ","p[i, D+1", aggregate_index, "] <- 1 - S[i, D", aggregate_index, "]\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:D) {\n")
      # Intialise linear predictor for delay
      linear_predictor_delay<-paste0(aggregate_space, "  ", "    ", link_function,"(S[i, d", aggregate_index, "]) <- alpha_0[d", aggregate_index, "] + ")

    }else{
      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
      # Calculate relative proportions nu from cumulative proportions S
      nimble_code <- paste0(nimble_code, aggregate_space, "  ","nu[i, 1", aggregate_index, "] <- S[i, 1", aggregate_index, "]\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 2:D) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "    ","nu[i, d", aggregate_index, "] <- (S[i, d", aggregate_index, "] - S[i, d-1", aggregate_index, "])/(1 - S[i, d-1", aggregate_index, "])\n")
      # End delay loop
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:D) {\n")
      # Intialise linear predictor for delay
      linear_predictor_delay<-paste0(aggregate_space, "  ", "    ", link_function,"(S[i, d", aggregate_index, "]) <- alpha_0[d", aggregate_index, "] + ")

    }

  } else if (stringr::str_detect(tolower(delay_link), "hazard")) {
    link_function <- stringr::str_extract(delay_link, '\\b\\w+$') # extract variable after "-"
    if(link_function=="hazard"){
      link_function<-"logit"
    }
    if(!(tolower(family$delay)=="gdm")){

      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
      # Calculate absolute proportions p from relative proportions nu
      nimble_code <- paste0(nimble_code, aggregate_space, "  ","p[i, 1", aggregate_index, "] <- nu[i, 1", aggregate_index, "]\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 2:D) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "    ","p[i, d", aggregate_index, "] <- (nu[i, d", aggregate_index, "])*(1 - sum(p[i, 1:(d-1)", aggregate_index, "])) \n")
      # End delay loop
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ","p[i, D+1", aggregate_index, "] <- (1- sum(p[i, 1:D", aggregate_index, "])) \n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:D) {\n")
      # Intialise linear predictor for delay
      linear_predictor_delay<-paste0(aggregate_space, "  ", "    ", link_function,"(nu[i, d", aggregate_index, "]) <- alpha_0[d", aggregate_index, "] + ")



    }else{
      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:D) {\n")
      # Intialise linear predictor for delay
      linear_predictor_delay<-paste0(aggregate_space, "  ","    ", link_function, "(nu[i, d", aggregate_index, "]) <- alpha_0[d", aggregate_index, "] + ")


    }

  } else if (stringr::str_detect(tolower(delay_link), "clr")) {
    if(D_index<2){
      stop(paste0("Please increase the argument for D such that it is greater than 2 or choose an alternative link function to 'clr'."))
    }
    link_function<-"clr"
    if(!(tolower(delay_link)=="clr")){
      print(paste0("Model version delay_link='clr' is being fit. Any additional information in delay_link is not recognised."))
    }
    if(!(tolower(family$delay)=="gdm")){
      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "p[i, D+1", aggregate_index, "] <- - sum(p[i, 1:D", aggregate_index, "])\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:D) {\n")
      # Model absolute proportions p
      # nimble_code <- paste0(nimble_code, aggregate_space, "  ","    ","nu[i, d", aggregate_index, "] <- exp(p[i, d", aggregate_index, "])/sum(exp(p[i, d:(D+1)", aggregate_index, "]))\n")
      linear_predictor_delay<-paste0(aggregate_space, "  ","    ","p[i, d", aggregate_index, "] <- alpha_0[d", aggregate_index, "] + ")

    }else{
      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "p[i, D+1", aggregate_index, "] <- - sum(p[i, 1:D", aggregate_index, "])\n")
      # Loop for delay
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:D) {\n")
      # Calculate relative proportions nu from absolute proportions p
      nimble_code <- paste0(nimble_code, aggregate_space, "  ","    ","nu[i, d", aggregate_index, "] <- exp(p[i, d", aggregate_index, "])/sum(exp(p[i, d:(D+1)", aggregate_index, "]))\n")
      linear_predictor_delay<-paste0(aggregate_space, "  ","    ","p[i, d", aggregate_index, "] <- alpha_0[d", aggregate_index, "] + ")

    }
  }else if (tolower(delay_link)=="logit"){
    link_function<-"logit"
    if(!(tolower(family$delay)=="nb"|tolower(family$delay)=="poisson")){
      stop("Argument delay_link='logit' is only valid if family$delay is 'nb' or 'poisson'.")
    }
    # Loop for time
    nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
    #nimble_code <- paste0(nimble_code, aggregate_space, "  ", "p[i, D+1", aggregate_index, "] <- 1 - sum(p[i, 1:D", aggregate_index, "])\n")
    # Loop for delay
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:(D+1)) {\n")
    # Calculate relative proportions nu from absolute proportions p
    linear_predictor_delay<-paste0(aggregate_space, "  ","    ","logit(p[i, d", aggregate_index, "]) <- alpha_0[d", aggregate_index, "] + ")


  }else if (tolower(delay_link)=="log"){
    link_function<-"log"
    if(!(family$delay=="nb"|tolower(family$delay)=="poisson")){
      stop("Argument delay_link='log' is only valid if family$delay is 'nb' or 'poisson'.")
    }
    # Loop for time
    nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
    # Loop for delay
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(d in 1:(D+1)) {\n")
    # Calculate relative proportions nu from absolute proportions p
    linear_predictor_delay<-paste0(aggregate_space, "  ","    ","log(mu[i, d", aggregate_index, "]) <- log(lambda[i",aggregate_index,"]) + alpha_0[d", aggregate_index, "] + ")

  }else{
    stop("Argument for delay_link is not recognised.")
  }

  if(!link_function%in%c("clr","log","logit","probit","cloglog")){
    print(paste0("Please note that the delay link function choice ",link_function ," may not be recommended/recognised."))
  }


  # List to store terms like splines, random effects, and structured terms
  spline_terms_delay <- NULL # to store spline term basis functions
  random_effects_delay <- list()  # for random effects
  structured_terms_delay <- list()  # for structured additive models
  jagam_output_delay  <- list()  # to store spline mgcv::jagam output

  jagam_output <- NULL
  for (pred_p in predictors_partial) {
    if(is.na(pred_p)){
      print(paste0("Model with just an intercept term in delay is being fit by default."))
    }else if (stringr::str_detect(pred_p, "s\\(")) {
      # Handle spline term with mgcv::jagam
      pred_name <- stringr::str_extract(pred_p, "(?<=\\().*(?=\\,)")  # extract variable inside s()

      # Use mgcv::jagam to get the spline basis for the variable
      data_spline_delay<-data
      data_spline_delay[response_partial]<-NULL
      data_spline_delay[[response_partial]]<-stats::rnorm(nimble_constants$N,0,1)
      jagam_output <- mgcv::jagam(stats::as.formula(paste(response_partial,"~",pred_p)), data = data_spline_delay, file='blank.jags')

      # Add the spline term in the linear pred_pictor
      spline_term_name <- paste0("g_spline_", pred_name)

      if(stringr::str_detect(tolower(delay_link), "survivor")){
        # If survivor link function, have same spline for each delay
        linear_predictor_delay <- paste0(linear_predictor_delay," alpha_spline_", pred_name, "[i", aggregate_index, "] + ", sep='')

        # Add spline basis to constants.
        nimble_constants[[paste("K_alpha_",pred_name,sep='')]]<-dim(jagam_output$jags.data$S1)[1]
        nimble_constants[[paste("S_alpha_",pred_name,sep='')]]<-jagam_output$jags.data$S1

        # Extract the spline basis matrix.
        # Calculate using mgcv::gam incase there is missing values:
        gam_output <- mgcv::gam(stats::as.formula(paste(response_partial,"~",pred_p)), data = data_spline_delay, fit=TRUE)
        blank_gam_x<-stats::predict(gam_output, newdata=data_spline_delay, type="lpmatrix")
        spline_basis <- blank_gam_x[,2:(dim(jagam_output$jags.data$S1)[1]+1)]
        #spline_basis <- jagam_output$jags.data$X[,2:(dim(jagam_output$jags.data$S1)[1]+1)]
        nimble_constants[[spline_term_name]]<-spline_basis

        # Add predictor to constants and define spline dimensions
        #nimble_constants[[pred_name]]<-data[[pred_name]]
        # if(aggregate){
        #   nimble_dimensions[[paste0("beta_spline_", pred_name)]]<-c(max(data[[pred_name]]),nimble_constants$A)
        #
        # }else(
        #   nimble_dimensions[[paste0("beta_spline_", pred_name)]]<-c(max(data[[pred_name]]))
        #
        # )

        # Store spline basis for NIMBLE
        spline_terms_delay[[spline_term_name]] <- pred_name
        jagam_output_delay[[spline_term_name]] <- jagam_output

        # Add initial values for the spline coefficients
        for(c in 1:nchains){
          if(aggregate){
            nimble_initial_values[[c]][[paste0("kappa_alpha_", pred_name)]] <- array(stats::rnorm(dim(spline_basis)[2]*nimble_constants$A, 0, 0.1), dim=c(dim(spline_basis)[2],nimble_constants$A))
          }else{
            nimble_initial_values[[c]][[paste0("kappa_alpha_", pred_name)]] <- stats::rnorm(dim(spline_basis)[2], 0, 0.1)
          }
        }

      }else{
        linear_predictor_delay <- paste0(linear_predictor_delay," alpha_spline_", pred_name, "[i, d", aggregate_index, "] + ", sep='')

        # Add spline basis to constants.
        nimble_constants[[paste("K_alpha_",pred_name,sep='')]]<-dim(jagam_output$jags.data$S1)[1]
        nimble_constants[[paste("S_alpha_",pred_name,sep='')]]<-jagam_output$jags.data$S1

        # Extract the spline basis matrix.
        # Calculate using mgcv::gam incase there is missing values:
        gam_output <- mgcv::gam(stats::as.formula(paste(response_partial,"~",pred_p)), data = data_spline_delay, fit=TRUE)
        blank_gam_x<-stats::predict(gam_output, newdata=data_spline_delay, type="lpmatrix")
        spline_basis <- blank_gam_x[,2:(dim(jagam_output$jags.data$S1)[1]+1)]
        #spline_basis <- jagam_output$jags.data$X[,2:(dim(jagam_output$jags.data$S1)[1]+1)]
        nimble_constants[[spline_term_name]]<-spline_basis

        # Add predictor to constants and define spline dimensions
        #nimble_constants[[pred_name]]<-data[[pred_name]]
        # if(aggregate){
        #   nimble_dimensions[[paste0("beta_spline_", pred_name)]]<-c(max(data[[pred_name]]),nimble_constants$A)
        #
        # }else(
        #   nimble_dimensions[[paste0("beta_spline_", pred_name)]]<-c(max(data[[pred_name]]))
        #
        # )

        # Store spline basis for NIMBLE
        spline_terms_delay[[spline_term_name]] <- pred_name
        jagam_output_delay[[spline_term_name]] <- jagam_output

        # Add initial values for the spline coefficients
        for(c in 1:nchains){
          if(aggregate){
            nimble_initial_values[[c]][[paste0("kappa_alpha_", pred_name)]] <- array(stats::rnorm(dim(spline_basis)[2]*nimble_constants$D*nimble_constants$A, 0, 0.1), dim=c(dim(spline_basis)[2], nimble_constants$D,nimble_constants$A))
          }else{
            nimble_initial_values[[c]][[paste0("kappa_alpha_", pred_name)]] <- matrix(stats::rnorm(dim(spline_basis)[2]*nimble_constants$D, 0, 0.1), ncol=nimble_constants$D)
          }
        }
      }

      nimble_monitors[[paste0("kappa_alpha_", pred_name)]] <- paste0("kappa_alpha_", pred_name)
      nimble_monitors[[paste0("alpha_spline_", pred_name)]] <- paste0("alpha_spline_", pred_name)


    } else if (stringr::str_detect(pred_p, "f\\(")) {
      # Handle random effects or structured models
      pred_name <- stringr::str_extract(pred_p, "(?<=\\().*(?=\\,)")  # extract variable inside f()
      if(stringr::str_detect(tolower(delay_link), "survivor")){
        # Determine the type of random effect or structured model
        if (stringr::str_detect(pred_p, "model = \"iid\"")) {
          # Random intercept
          nimble_constants[[paste0(pred_name,"_factor")]]<-as.numeric(as.factor(data[[pred_name]]))
          nimble_constants[[paste0(pred_name,"_length")]]<-length(unique(nimble_constants[[paste0(pred_name,"_factor")]]))
          linear_predictor_delay <- paste0(linear_predictor_delay, "v_", pred_name, "[",pred_name,"_factor[i]", aggregate_index, "] + ")
          random_effects_delay[[pred_name]] <- "iid"

          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste0("tav_", pred_name)]] <- abs(stats::rnorm(nimble_constants$A,0,1))
              nimble_initial_values[[c]][[paste0("v_", pred_name)]] <- array((stats::rnorm(nimble_constants$A*nimble_constants$N,0,1)), dim=c(nimble_constants$N,nimble_constants$A))

            }else{
              nimble_initial_values[[c]][[paste0("tav_", pred_name)]] <- abs(stats::rnorm(1,0,1))
              nimble_initial_values[[c]][[paste0("v_", pred_name)]] <- stats::rnorm(nimble_constants$N,0,1)

            }
          }
          nimble_monitors[[paste0("tav_", pred_name)]] <- paste0("tav_", pred_name)
          nimble_monitors[[paste0("v_", pred_name)]] <- paste0("v_", pred_name)

        } else if (stringr::str_detect(pred_p, "model = \"rw1\"")) {
          # Random walk (RW1)
          linear_predictor_delay <- paste0(linear_predictor_delay, "v_rw1_", pred_name, "[i", aggregate_index, "] + ")
          structured_terms_delay[[pred_name]] <- "rw1"
          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste0("tav_rw1_", pred_name)]] <- abs(stats::rnorm(nimble_constants$A,0,1))
              nimble_initial_values[[c]][[paste0("v_rw1_", pred_name)]] <- array((stats::rnorm(nimble_constants$A*nimble_constants$N,0,1)), dim=c(nimble_constants$N,nimble_constants$A))

            }else{
              nimble_initial_values[[c]][[paste0("tav_rw1_", pred_name)]] <- abs(stats::rnorm(1,0,1))
              nimble_initial_values[[c]][[paste0("v_rw1_", pred_name)]] <- stats::rnorm(nimble_constants$N,0,1)

            }
          }
          nimble_monitors[[paste0("tav_rw1_", pred_name)]] <- paste0("tav_rw1_", pred_name)
          nimble_monitors[[paste0("v_rw1_", pred_name)]] <- paste0("v_rw1_", pred_name)

        }
      }else{
        if (stringr::str_detect(pred_p, "model = \"iid\"")) {
          # Random intercept
          nimble_constants[[paste0(pred_name,"_factor")]]<-as.numeric(as.factor(data[[pred_name]]))
          nimble_constants[[paste0(pred_name,"_length")]]<-length(unique(nimble_constants[[paste0(pred_name,"_factor")]]))
          linear_predictor_delay <- paste0(linear_predictor_delay, "v_", pred_name, "[",pred_name,"_factor[i], d", aggregate_index, "] + ")
          random_effects_delay[[pred_name]] <- "iid"

          for(c in 1:nchains){
            if ((tolower(delay_link)=="log"|(tolower(delay_link)=="logit"))){
              if(aggregate){
                nimble_initial_values[[c]][[paste0("tav_", pred_name)]] <- matrix(abs(stats::rnorm(nimble_constants$A*(nimble_constants$D+1),0,1)), ncol=nimble_constants$A)
                nimble_initial_values[[c]][[paste0("v_", pred_name)]] <- array((stats::rnorm(nimble_constants$A*(nimble_constants$D+1)*nimble_constants$N,0,1)), dim=c(nimble_constants$N,(nimble_constants$D+1),nimble_constants$A))

              }else{
                nimble_initial_values[[c]][[paste0("tav_", pred_name)]] <- abs(stats::rnorm((nimble_constants$D+1),0,1))
                nimble_initial_values[[c]][[paste0("v_", pred_name)]] <- matrix((stats::rnorm((nimble_constants$D+1)*nimble_constants$N,0,1)), ncol=(nimble_constants$D+1))

              }
            }else{
              if(aggregate){
                nimble_initial_values[[c]][[paste0("tav_", pred_name)]] <- matrix(abs(stats::rnorm(nimble_constants$A*nimble_constants$D,0,1)), ncol=nimble_constants$A)
                nimble_initial_values[[c]][[paste0("v_", pred_name)]] <- array((stats::rnorm(nimble_constants$A*nimble_constants$D*nimble_constants$N,0,1)), dim=c(nimble_constants$N,nimble_constants$D,nimble_constants$A))

              }else{
                nimble_initial_values[[c]][[paste0("tav_", pred_name)]] <- abs(stats::rnorm(nimble_constants$D,0,1))
                nimble_initial_values[[c]][[paste0("v_", pred_name)]] <- matrix((stats::rnorm(nimble_constants$D*nimble_constants$N,0,1)), ncol=nimble_constants$D)

              }            }

          }
          nimble_monitors[[paste0("tav_", pred_name)]] <- paste0("tav_", pred_name)
          nimble_monitors[[paste0("v_", pred_name)]] <- paste0("v_", pred_name)


        } else if (stringr::str_detect(pred_p, "model = \"rw1\"")) {
          # Random walk (RW1)
          linear_predictor_delay <- paste0(linear_predictor_delay, "v_rw1_", pred_name, "[i, d", aggregate_index, "] + ")
          structured_terms_delay[[pred_name]] <- "rw1"
          for(c in 1:nchains){
            if ((tolower(delay_link)=="log"|(tolower(delay_link)=="logit"))){
              if(aggregate){
                nimble_initial_values[[c]][[paste0("tav_rw1_", pred_name)]] <- matrix(abs(stats::rnorm(nimble_constants$A*(nimble_constants$D+1),0,1)), ncol=nimble_constants$A)
                nimble_initial_values[[c]][[paste0("v_rw1_", pred_name)]] <- array((stats::rnorm(nimble_constants$A*(nimble_constants$D+1)*nimble_constants$N,0,1)), dim=c(nimble_constants$N,(nimble_constants$D+1),nimble_constants$A))

              }else{
                nimble_initial_values[[c]][[paste0("tav_rw1_", pred_name)]] <- abs(stats::rnorm((nimble_constants$D+1),0,1))
                nimble_initial_values[[c]][[paste0("v_rw1_", pred_name)]] <- matrix((stats::rnorm((nimble_constants$D+1)*nimble_constants$N,0,1)), ncol=(nimble_constants$D+1))

              }
            }else{
              if(aggregate){
                nimble_initial_values[[c]][[paste0("tav_rw1_", pred_name)]] <- matrix(abs(stats::rnorm(nimble_constants$A*nimble_constants$D,0,1)), ncol=nimble_constants$A)
                nimble_initial_values[[c]][[paste0("v_rw1_", pred_name)]] <- array((stats::rnorm(nimble_constants$A*nimble_constants$D*nimble_constants$N,0,1)), dim=c(nimble_constants$N,nimble_constants$D,nimble_constants$A))

              }else{
                nimble_initial_values[[c]][[paste0("tav_rw1_", pred_name)]] <- abs(stats::rnorm(nimble_constants$D,0,1))
                nimble_initial_values[[c]][[paste0("v_rw1_", pred_name)]] <- matrix((stats::rnorm(nimble_constants$D*nimble_constants$N,0,1)), ncol=nimble_constants$D)

              }
            }

          }
          nimble_monitors[[paste0("tav_rw1_", pred_name)]] <- paste0("tav_rw1_", pred_name)
          nimble_monitors[[paste0("v_rw1_", pred_name)]] <- paste0("v_rw1_", pred_name)

        }
      }

    } else if (pred_p == "1"|pred_p == "I") {
      # Intercept included by default.
      # if(stringr::str_detect(tolower(delay_link), "survivor")){
      #   # Intercept term
      #   linear_predictor_delay <- paste0(linear_predictor_delay, "alpha_int", aggregate_single," + ")
      #   for(c in 1:nchains){
      #   if(aggregate){
      #     nimble_initial_values[[c]][["alpha_int"]] <- stats::rnorm(nimble_constants$A, 0, 1)
      #   }else{
      #     nimble_initial_values[[c]][["alpha_int"]] <- stats::rnorm(1, 0, 1)
      #   }
      #   }
      #
      # }else{
      #   # Intercept term
      #   linear_predictor_delay <- paste0(linear_predictor_delay, "alpha_int[d", aggregate_index, "] + ")
      #   for(c in 1:nchains){
      #   if(aggregate){
      #     nimble_initial_values[[c]][["alpha_int"]] <- matrix(stats::rnorm(nimble_constants$A*nimble_constants$D, 0, 1), ncol=nimble_constants$A)
      #   }else{
      #     nimble_initial_values[[c]][["alpha_int"]] <- stats::rnorm(nimble_constants$D, 0, 1)
      #   }
      #   }
      #
      # }
      # nimble_monitors[[paste0("alpha_int")]] <- "alpha_int"
    } else {
      # Add covariate to data
      nimble_data[[pred_p]]<-data[[pred_p]]
      if(is.null(data[[pred_p]])){
        stop(paste0("Please include covariate ", pred, " in data argument."))
      }
      if(stringr::str_detect(tolower(delay_link), "survivor")){
        # Regular predictor
        linear_predictor_delay <- paste0(linear_predictor_delay, "alpha_", pred_p, aggregate_single," * ", pred_p, "[i", aggregate_index, "] + ")
        for(c in 1:nchains){
          if(aggregate){
            nimble_initial_values[[c]][[paste0("alpha_", pred_p)]] <- stats::rnorm(nimble_constants$A, 0, 0.1)
          }else{
            nimble_initial_values[[c]][[paste0("alpha_", pred_p)]] <- stats::rnorm(1, 0, 0.1)
          }
        }
      }else{
        # Regular predictor
        linear_predictor_delay <- paste0(linear_predictor_delay, "alpha_", pred_p, "[d", aggregate_index, "] * ", pred_p, "[i", aggregate_index, "] + ")
        # EDIT above assumes delay covariate is just over time (is this more likely?)? - have separate time covariate for delay and time ?
        for(c in 1:nchains){
          if(aggregate){
            nimble_initial_values[[c]][[paste0("alpha_", pred_p)]] <- matrix(stats::rnorm(nimble_constants$A*nimble_constants$D, 0, 0.1), ncol=nimble_constants$A)
          }else{
            nimble_initial_values[[c]][[paste0("alpha_", pred_p)]] <- stats::rnorm(nimble_constants$D, 0, 0.1)
          }
        }
      }
      nimble_monitors[[paste0("alpha_", pred_p)]] <- paste0("alpha_", pred_p)

    }
  }

  # Remove trailing " + " from the linear predictor
  linear_predictor_delay <- substr(linear_predictor_delay, 1, nchar(linear_predictor_delay) - 3)

  # Add the linear predictor to the code
  nimble_code <- paste0(nimble_code, linear_predictor_delay, "\n")
  nimble_code <- paste0(nimble_code, aggregate_space,  "  ", "}\n")
  nimble_code <- paste0(nimble_code, aggregate_space, "}\n")

  ############ Model for nested counts:

  if(nested){

    # Parse the formula$nested components
    response_nested <- as.character(formula$nested[[2]])  # Dependent variable(s)
    predictors_nested <- as.character(formula$nested[-2])  # Independent variables

    # Clean up predictors (removes "+" and spaces)
    predictors_nested <- stringr::str_split(predictors_nested[2], "\\+")[[1]]
    predictors_nested <- stringr::str_trim(predictors_nested)

    # Add the likelihood model
    nimble_monitors[[paste0(response_nested, "_corrected")]] <- paste0(response_nested, "_corrected")
    nimble_monitors[[paste0("chi")]] <- paste0("chi")
    nimble_monitors[[paste0("eta")]] <- paste0("eta")
    for(c in 1:nchains){
      if(aggregate){
        # EDIT FOR multiple virus: nimble_initial_values[[c]][[paste0("chi")]] <- matrix(abs(stats::rnorm(nimble_constants$A, 0, 5)), ncol=nimble_constants$A)
        nimble_initial_values[[c]][[paste0("chi")]] <- abs(stats::rnorm(nimble_constants$A, 0, 5))
      }else{
        nimble_initial_values[[c]][[paste0("chi")]] <- abs(stats::rnorm(1, 0, 5))
      }
    }
    # Loop for time
    nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:W) {\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", response_nested, "_corrected[i", aggregate_index, "] ~ ")
    nimble_code <- paste0(nimble_code, "dbetabin(eta[i", aggregate_index, "], chi", aggregate_single, ", ", response, "[i", aggregate_index, "])\n")
    # End time loop
    nimble_code <- paste0(nimble_code, aggregate_space,  "}\n")







    # Create linear predictor for response for the partial counts

    # EDIT - if(nested.multiple){} - multiple virus than need same as GDM link options
    # } else if(tolower(link$nested) == "survivor") {
    #   # Loop for time
    #   nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
    #   # Calculate relative proportions eta from cumulative proportions S
    #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "eta[i", aggregate_index, "] <- (S[i", aggregate_index, "] - S[i", aggregate_index, "])/(1 - S[i", aggregate_index, "])\n")
    #   # Intialise linear predictor for nested
    #   linear_predictor_nested<-paste0(aggregate_space, "  ", "probit(S[i", aggregate_index, "]) <- gamma_0", aggregate_single, " +")
    #
    # } else if (tolower(link$nested) == "hazard") {
    #   # Loop for time
    #   nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
    #   # Intialise linear predictor for nested
    #   linear_predictor_nested<-paste0(aggregate_space, "  ", "logit(eta[i", aggregate_index, "]) <- gamma_0",aggregate_single," +")
    #
    # } else if (tolower(link$nested) == "clr") {
    #   # Loop for time
    #   nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
    #   # Calculate relative proportions eta from absolute proportions p
    #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "eta[i", aggregate_index, "] <- exp(p[i", aggregate_index, "])/sum(exp(p[i, d:(D+1)", aggregate_index, "]))\n")
    #   linear_predictor_nested<-paste0(aggregate_space, "  " ,"p[i, d", aggregate_index, "] <- gamma_0[d", aggregate_index, "] +")
    #

    # Loop for time
    nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:N) {\n")
    # Initialise linear predictor for nested
    linear_predictor_nested<-paste0(aggregate_space, "  ", "logit(eta[i", aggregate_index, "]) <- gamma_0",aggregate_single," + delta",aggregate_single, "*( log(lambda[i", aggregate_index, "]) - beta_0", aggregate_single, ") + ")


    # List to store terms like splines, random effects, and structured terms
    spline_terms_nested <- NULL # to store spline term basis functions
    random_effects_nested <- list()  # for random effects
    structured_terms_nested <- list()  # for structured additive models
    jagam_output_nested  <- list()  # to store spline mgcv::jagam output

    jagam_output <- NULL
    for (pred_p in predictors_nested) {
      if (stringr::str_detect(pred_p, "s\\(")) {
        # Handle spline term with mgcv::jagam
        pred_name <- stringr::str_extract(pred_p, "(?<=\\().*(?=\\,)")  # extract variable inside s()

        # Use mgcv::jagam to get the spline basis for the variable
        data_spline_nested<-data
        data_spline_nested[response_nested]<-NULL
        data_spline_nested[[response_nested]]<-stats::rnorm(nimble_constants$N,0,1)
        jagam_output <- mgcv::jagam(stats::as.formula(paste(response_nested,"~",pred_p)), data = data_spline_nested, file='blank.jags')

        # Add the spline term in the linear pred_pictor
        spline_term_name <- paste0("g_spline_", pred_name)
        linear_predictor_nested <- paste0(linear_predictor_nested," gamma_spline_", pred_name, "[i", aggregate_index, "] + ", sep='')

        # Add spline basis to constants.
        nimble_constants[[paste("K_gamma_",pred_name,sep='')]]<-dim(jagam_output$jags.data$S1)[1]
        nimble_constants[[paste("S_gamma_",pred_name,sep='')]]<-jagam_output$jags.data$S1

        # Extract the spline basis matrix.
        # spline_basis <- jagam_output$jags.data$X[,2:(dim(jagam_output$jags.data$S1)[1]+1)]
        gam_output <- mgcv::gam(stats::as.formula(paste(response_nested,"~",pred_p)), data = data_spline_nested, fit=TRUE)
        blank_gam_x<-stats::predict(gam_output, newdata=data_spline_nested, type="lpmatrix")
        spline_basis <- blank_gam_x[,2:(dim(jagam_output$jags.data$S1)[1]+1)]

        nimble_constants[[spline_term_name]]<-spline_basis

        # Add predictor to constants and define spline dimensions
        #nimble_constants[[pred_name]]<-data[[pred_name]]
        # if(aggregate){
        #   nimble_dimensions[[paste0("beta_spline_", pred_name)]]<-c(max(data[[pred_name]]),nimble_constants$A)
        #
        # }else(
        #   nimble_dimensions[[paste0("beta_spline_", pred_name)]]<-c(max(data[[pred_name]]))
        #
        # )

        # Store spline basis for NIMBLE
        spline_terms_nested[[spline_term_name]] <- pred_name
        jagam_output_nested[[spline_term_name]] <- jagam_output

        # Add initial values for the spline coefficients
        for(c in 1:nchains){
          if(aggregate){
            nimble_initial_values[[c]][[paste0("kappa_gamma_", pred_name)]] <- matrix(stats::rnorm(dim(spline_basis)[2]*nimble_constants$A, 0, 0.1), ncol=nimble_constants$A)
          }else{
            nimble_initial_values[[c]][[paste0("kappa_gamma_", pred_name)]] <- stats::rnorm(dim(spline_basis)[2], 0, 0.1)
          }
        }

        nimble_monitors[[paste0("kappa_gamma_", pred_name)]] <- paste0("kappa_gamma_", pred_name)
        nimble_monitors[[paste0("gamma_spline_", pred_name)]] <- paste0("gamma_spline_", pred_name)

      } else if (stringr::str_detect(pred_p, "f\\(")) {
        # Handle random effects or structured models
        pred_name <- stringr::str_extract(pred_p, "(?<=\\().*(?=\\,)")  # extract variable inside f()

        # Determine the type of random effect or structured model
        if (stringr::str_detect(pred_p, "model = \"iid\"")) {
          # Random intercept
          nimble_constants[[paste0(pred_name,"_factor")]]<-as.numeric(as.factor(data[[pred_name]]))
          nimble_constants[[paste0(pred_name,"_length")]]<-length(unique(nimble_constants[[paste0(pred_name,"_factor")]]))
          linear_predictor_nested <- paste0(linear_predictor_nested, "w_", pred_name, "[",pred_name,"_factor[i]", aggregate_index, "] + ")
          random_effects_nested[[pred_name]] <- "iid"
          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste0("taw_", pred_name)]] <- abs(stats::rnorm(nimble_constants$A,0,1))
            }else{
              nimble_initial_values[[c]][[paste0("taw_", pred_name)]] <- abs(stats::rnorm(1,0,1))
            }
          }
          nimble_monitors[[paste0("w_", pred_name)]] <- paste0("w_", pred_name)
          nimble_monitors[[paste0("taw_", pred_name)]] <- paste0("taw_", pred_name)

        } else if (stringr::str_detect(pred_p, "model = \"rw1\"")) {
          # Random walk (RW1)
          linear_predictor_nested <- paste0(linear_predictor_nested, "w_rw1_", pred_name, "[i", aggregate_index, "] + ")
          structured_terms_nested[[pred_name]] <- "rw1"
          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste0("taw_rw1_", pred_name)]] <- abs(stats::rnorm(nimble_constants$A,0,1))
            }else{
              nimble_initial_values[[c]][[paste0("taw_rw1_", pred_name)]] <- abs(stats::rnorm(1,0,1))
            }
          }
          nimble_monitors[[paste0("w_rw1_", pred_name)]] <- paste0("w_rw1_", pred_name)
          nimble_monitors[[paste0("taw_rw1_", pred_name)]] <- paste0("taw_rw1_", pred_name)

          # } else if (stringr::str_detect(pred_p, "model = \"rw2\"")) {
          #   # Random walk 2 (RW2)
          #   linear_predictor_nested <- paste0(linear_predictor_nested, "w_rw2_", pred_name, "[i", aggregate_index, "] + ")
          #   structured_terms_nested[[pred_name]] <- "rw2"
          #   for(c in 1:nchains){
          #   if(aggregate){
          #     nimble_initial_values[[c]][[paste0("taw_rw2_", pred_name)]] <-  abs(stats::rnorm(nimble_constants$A,0,1))
          #   }else{
          #     nimble_initial_values[[c]][[paste0("taw_rw2_", pred_name)]] <-  abs(stats::rnorm(1,0,1))
          #   }
          #   }
          #   nimble_monitors[[paste0("w_rw2_", pred_name)]] <- paste0("w_rw2_", pred_name)
          #   nimble_monitors[[paste0("taw_rw2_", pred_name)]] <- paste0("taw_rw2_", pred_name)
          #
          #   # } else if (stringr::str_detect(pred_p, "model = \"seasonal\"")) {
          #   #   # Seasonality (periodic effect)
          #   #   linear_predictor_nested <- paste0(linear_predictor_nested, "u_seasonal_", pred_name, "[i", aggregate_index, "] + ")
          # } else if (stringr::str_detect(pred_p, "model = \"spatial\"")) {
          #   # Spatial effect (spatial random field)
          #   linear_predictor_nested <- paste0(linear_predictor_nested, "w_spatial[i", aggregate_index, "] + ")
          #   structured_terms_nested[["spatial"]] <- "spatial"
          #   for(c in 1:nchains){
          #   if(aggregate){
          #     nimble_initial_values[[c]][["taw_spatial"]] <- abs(stats::rnorm(nimble_constants$A,0,1))
          #   }else{
          #     nimble_initial_values[[c]][["taw_spatial"]] <- abs(stats::rnorm(1,0,1))
          #   }
          #   }
          # }
          # nimble_monitors[[paste0("taw_spatial")]] <- paste0("taw_spatial")
          # nimble_monitors[[paste0("w_spatial")]] <- paste0("w_spatial")
          #
        } else if (pred_p == "1"|pred_p == "I") {
          print("Please note an intercept term is included in the nested formula by default.")
        } else {
          # Regular predictor
          linear_predictor_nested <- paste0(linear_predictor_nested, "gamma_", pred_p, aggregate_single, " * ", pred_p, "[i", aggregate_index, "] + ")
          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste0("gamma_", pred_p)]] <- stats::rnorm(nimble_constants$A, 0, 0.1)
            }else{
              nimble_initial_values[[c]][[paste0("gamma_", pred_p)]] <- stats::rnorm(1, 0, 0.1)
            }
          }
          nimble_monitors[[paste0("gamma_", pred_p)]] <- paste0("gamma_", pred_p)

        }
      }

      # Remove trailing " + " from the linear predictor
      linear_predictor_nested <- substr(linear_predictor_nested, 1, nchar(linear_predictor_nested) - 3)

      # Add the linear predictor to the code
      nimble_code <- paste0(nimble_code, linear_predictor_nested, "\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "}\n")

      # Model for censored nested variable

      if(is.null(nested.censored)){
        print(paste0("Argument nested.censored has not been set so nested variable is assumed to be fully observed."))
        nested.censored <- 0
      }
      # Constant for period up to which nested variable is not censored.
      nimble_constants$C <- nimble_constants$N -  nested.censored

      # Add the likelihood model
      for(c in 1:nchains){
        if(aggregate){
          nimble_initial_values[[c]][[paste0("upsilon")]] <-abs(stats::rnorm(nimble_constants$A, 0, 5))
        }else{
          nimble_initial_values[[c]][[paste0("upsilon")]] <- abs(stats::rnorm(1, 0, 5))
        }
      }
      nimble_monitors[[paste0("upsilon")]] <- paste0("upsilon")

      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in (C+1):W) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", response_nested, "[i", aggregate_index, "] ~ ")
      nimble_code <- paste0(nimble_code, "dbetabin(pi[i", aggregate_index, "], upsilon", aggregate_single, ", ", response_nested, "_corrected[i", aggregate_index, "])\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "pi[i", aggregate_index, "] <- exp((-omega",aggregate_single,")*(i-C)) \n")
      # End time loop
      nimble_code <- paste0(nimble_code, aggregate_space,  "}\n")

      ## When nested proportion is not censored pi = 1.
      # Loop for time
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:C) {\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "pi[i", aggregate_index, "] <- 1 \n")
      # End time loop
      nimble_code <- paste0(nimble_code, aggregate_space, "}\n")

      nimble_monitors[[paste0("pi")]] <- paste0("pi")
      nimble_monitors[[paste0("omega")]] <- paste0("omega")
      nimble_monitors[[paste0("delta")]] <- paste0("delta")



    }
  }


  ############ Priors:
  ### Set prior distributions
  # Set defaults:
  priors_default<-list(
    # default priors for totals
    totals.intercept="dnorm(0, 10)",
    totals.dispersion="dgamma(2,0.02)",
    totals.spline="T(dnorm(0,1),0,)",
    totals.rw2="T(dnorm(0,1),0,)",
    totals.rw1="T(dnorm(0,1),0,)",
    totals.iid="T(dnorm(0,1),0,)",
    totals.re="dnorm(0, 10)",

    # default priors for partial counts
    delay.intercept="dnorm(0, 10)",
    delay.dispersion="dgamma(2,0.02) ",
    delay.spline="T(dnorm(0,1),0,)",
    delay.rw2="T(dnorm(0,1),0,)",
    delay.rw1="T(dnorm(0,1),0,)",
    delay.iid="T(dnorm(0,1),0,)",
    delay.re="dnorm(0, 10)",

    # default priors for nested counts
    nested.intercept="dnorm(0, 10)",
    nested.dispersion="dgamma(2,0.02)",
    nested.scale="dnorm(0,sd=5)",
    nested.spline="T(dnorm(0,1),0,)",
    nested.rw2="T(dnorm(0,1),0,)",
    nested.rw1="T(dnorm(0,1),0,)",
    nested.iid="T(dnorm(0,1),0,)",
    nested.re="dnorm(0, 10)",
    censored.dispersion="dgamma(2,0.02)",
    censored.slope="dgamma(shape=10,rate=200)"
  )
  nimble_priors<-priors_default
  # Replace any defaults with distributions set in priors argument:
  nimble_priors[names(priors)]<-priors



  ### Priors for totals model

  # Add prior distributions for the coefficients
  nimble_code <- paste0(nimble_code, aggregate_space, "# Priors total\n")
  nimble_code <- paste0(nimble_code, aggregate_space, "beta_0", aggregate_single, " ~ ",nimble_priors$totals.intercept,"\n")
  for(c in 1:nchains){
    if(aggregate){
      nimble_initial_values[[c]][["beta_0"]]<- stats::rnorm(nimble_constants$A,0,0.1)
    }else{
      nimble_initial_values[[c]][["beta_0"]]<- stats::rnorm(1,0,0.1)
    }
  }
  nimble_monitors[[paste0("beta_0")]]<-paste0("beta_0")

  # Add priors for regular predictors
  for (pred in predictors) {
    if (!stringr::str_detect(pred, "s\\(") && pred != "1" && !stringr::str_detect(pred, "f\\(")) {
      nimble_code <- paste0(nimble_code, aggregate_space, "beta_", pred, aggregate_single, " ~ ",nimble_priors$totals.re,"\n")

    }
  }

  if(!is.null(spline_terms)|!is.null(spline_terms_delay)){
    nimble_constants$zeros<-c(0)
  }
  # Add priors for spline coefficients
  if (!is.null(spline_terms)) {
    nimble_code <- paste0(nimble_code, aggregate_space, "# Priors for totals spline \n")
    for(j in 1:length(spline_terms)){
      # nimble_constants[[paste("S_beta_",pred_name,i,sep='')]]<-jagam_output$jags.data$S1[,(1+(i-1)*dim(jagam_output$jags.data$S1)[1]):(i*dim(jagam_output$jags.data$S1)[1])]
      nimble_code <- paste0(nimble_code, aggregate_space, "beta_spline_",spline_terms[j],"[1:N", aggregate_index, "] <- ","f_spline_",spline_terms[j], "[1:N, 1:K_beta_",pred_name,"]%*%kappa_beta_", pred_name, "[1:K_beta_",pred_name, aggregate_index, "] \n")
      nimble_code <- paste0(nimble_code, aggregate_space, "kappa_beta_",spline_terms[j],"[", "1:K_beta_",spline_terms[j], aggregate_index, "] ~ dmnorm(zeros[", "1:K_beta_",spline_terms[j],"],","omega_beta_",spline_terms[j],"[", "1:K_beta_",spline_terms[j], ", 1:K_beta_",spline_terms[j], aggregate_index, "])\n",sep='')
      nimble_constants$zeros<-rep(0,max(c(length(nimble_constants$zeros), nimble_constants[[paste("K_beta_",spline_terms[j],sep='')]])))
      # Penalty matrix for MVN distribution.
      nimble_code <- paste0(nimble_code, aggregate_space, "omega_beta_",spline_terms[j],"[", "1:K_beta_",spline_terms[j], ", 1:K_beta_",spline_terms[j], aggregate_index, "] <-",sep='')
      S_n<-(dim(jagam_output_totals[[j]]$jags.data$S1)[2]/dim(jagam_output_totals[[j]]$jags.data$S1)[1])
      for(i in 1:S_n){
        nimble_code <- paste0(nimble_code,"S_beta_",spline_terms[j],"[", "1:K_beta_",spline_terms[j], ", (1 + ",(i-1),"*K_beta_",spline_terms[j],"):(",i,"*K_beta_",spline_terms[j], ")]/sigma_beta_",spline_terms[j],i,aggregate_single,"^2 + ",sep='')
        if(i==S_n){
          nimble_code <- substr(nimble_code, 1, nchar(nimble_code) - 3)
          nimble_code <- paste0(nimble_code," \n",sep='')
        }
      }
      for(i in 1:S_n){
        nimble_code <- paste0(nimble_code, aggregate_space, "sigma_beta_", spline_terms[j],i,aggregate_single," ~ ",nimble_priors$totals.spline," \n")
        for(c in 1:nchains){
          if(aggregate){
            nimble_initial_values[[c]][[paste("sigma_beta_", spline_terms[j],i,sep='')]]<- stats::runif(nimble_constants$A, 0, 1)
          }else{
            nimble_initial_values[[c]][[paste("sigma_beta_", spline_terms[j],i,sep='')]]<- stats::runif(1, 0, 1)
          }
        }
        nimble_monitors[[paste("sigma_beta_", spline_terms[j],i,sep='')]]<-paste("sigma_beta_", spline_terms[j],i,sep='')

      }
    }
  }

  # Add priors for random effects and structured terms
  if (length(random_effects) > 0) {
    for (random_name in names(random_effects)) {
      nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:",random_name,"_length){\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "  ", "u_", random_name, "[i",aggregate_index,"] ~ dnorm(0, tau_", random_name, aggregate_single, ")\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
      nimble_code <- paste0(nimble_code, aggregate_space, "tau_", random_name, aggregate_single, " ~  ",nimble_priors$totals.iid, "\n")

    }
  }

  if (length(structured_terms) > 0) {
    for (term_name in names(structured_terms)) {
      if (structured_terms[[term_name]] == "rw1") {
        nimble_code <- paste0(nimble_code, aggregate_space, "u_rw1_", term_name, "[1", aggregate_index, "] ~ dnorm(0",", ", "tau_rw1_", term_name, aggregate_single, ")\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "for(j_", term_name, " in 2:N){\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "  ", "u_rw1_", term_name, "[j_", term_name, aggregate_index, "] ~ dnorm(u_rw1_", term_name,"[j_", term_name, "-1", aggregate_index,"]",", ", "tau_rw1_", term_name, aggregate_single, ")\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "tau_rw1_", term_name, aggregate_single, " ~ ", nimble_priors$totals.rw1, "\n")
      } else if (structured_terms[[term_name]] == "rw2") {
        nimble_code <- paste0(nimble_code, aggregate_space, "u_rw2_", term_name, "[1", aggregate_index,"] ~ dnorm(0, tau_rw2_", term_name, aggregate_single, ")\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "u_rw2_", term_name, "[2", aggregate_index,"] ~ dnorm(u_rw2_", term_name, "[1",aggregate_index,"], ", "tau_rw2_", term_name, aggregate_single, ")\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "for(j_", term_name, " in 3:N){\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "    ", "u_rw2_", term_name, "[j_", term_name, aggregate_index , "] ~ dnorm(2*u_rw2_", term_name,"[j_", term_name, "-1", aggregate_index,"] - u_rw2_", term_name, "[j_", term_name, "-2", aggregate_index,"]", ", ", "tau_rw2_", term_name, aggregate_single, ")\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "tau_rw2_", term_name, aggregate_single, " ~ ", nimble_priors$totals.rw2, "\n")
        # } else if (structured_terms[[term_name]] == "seasonal") {
        #   nimble_code <- paste0(nimble_code, "  u_seasonal_", term_name, "[1:N] ~ dPeriodic(tau_seasonal_", term_name, ")\n")
        #   nimble_code <- paste0(nimble_code, "  tau_seasonal_", term_name, " ~ T(dnorm(0,1),0,)\n")
      }
      # else if (structured_terms[[term_name]] == "spatial") { # EDIT allow for spatial effect
      #   nimble_code <- paste0(nimble_code, aggregate_space, "u_spatial[1:N", aggregate_index,"] ~ dcar_normal(adj[], num[], tau_spatial",aggregate_single,")\n")
      #   nimble_code <- paste0(nimble_code, aggregate_space, "tau_spatial",aggregate_single," ~ T(dnorm(0,1),0,)\n")
      # }
    }
  }

  # Add priors for family-specific parameters
  if(!is.null(family$total)){
    if ((tolower(family$total)=="nb")){
      nimble_code <- paste0(nimble_code, aggregate_space, "theta", aggregate_single, " ~ ", nimble_priors$totals.dispersion, "\n")
    }
  }

  ### Priors for delay model

  # Add prior distributions for the coefficients
  nimble_code <- paste0(nimble_code, aggregate_space, "# Priors delay\n")
  if (stringr::str_detect(tolower(delay_link), "survivor")) {
    nimble_code <- paste0(nimble_code, aggregate_space, "alpha_0[1", aggregate_index,"] ~ ", nimble_priors$delay.intercept, "\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 2:D){\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  alpha_0[d", aggregate_index,"] ~ T(dnorm(alpha_0[d-1", aggregate_index,"], sd=10), alpha_0[d-1", aggregate_index,"], )\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][["alpha_0"]]<- matrix(sort(stats::rnorm(nimble_constants$D*nimble_constants$A,0,0.1)), ncol=nimble_constants$A, byrow = TRUE)
      }else{
        nimble_initial_values[[c]][["alpha_0"]]<- sort(stats::rnorm(nimble_constants$D,0,0.1))
      }
    }


  } else if ((tolower(delay_link)=="log")|(tolower(delay_link)=="logit")){ # EDIT - also multinomial?
    nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:(D+1)){\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "alpha_0[d", aggregate_index,"] ~ ", nimble_priors$delay.intercept, "\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][["alpha_0"]]<- matrix(sort(stats::rnorm((nimble_constants$D+1)*nimble_constants$A,0,0.1)), ncol=nimble_constants$A, byrow = TRUE)
      }else{
        nimble_initial_values[[c]][["alpha_0"]]<- sort(stats::rnorm((nimble_constants$D+1),0,0.1))
      }
    }


  }else{

    nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:D){\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  alpha_0[d", aggregate_index,"] ~ ", nimble_priors$delay.intercept, "\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][["alpha_0"]]<- matrix(stats::rnorm(nimble_constants$D*nimble_constants$A,0,0.1), ncol=nimble_constants$A)
      }else{
        nimble_initial_values[[c]][["alpha_0"]]<- stats::rnorm(nimble_constants$D,0,0.1)
      }
    }

  }


  nimble_monitors[[paste0("alpha_0")]]<-paste0("alpha_0")


  # Add priors for regular predictors_partial
  for (pred_p in predictors_partial) {
    if(stringr::str_detect(tolower(delay_link), "survivor")){
      if (!stringr::str_detect(pred_p, "s\\(") && pred_p != "1" && !stringr::str_detect(pred_p, "f\\(")) {
        nimble_code <- paste0(nimble_code, aggregate_space, "  ", "alpha_", pred_p, aggregate_single," ~ ", nimble_priors$delay.re, "\n")
        for(c in 1:nchains){
          if(aggregate){
            nimble_initial_values[[c]][[paste0("alpha_", pred_p)]]<- stats::rnorm(nimble_constants$A,0,0.1)
          }else{
            nimble_initial_values[[c]][[paste0("alpha_", pred_p)]]<- stats::rnorm(1,0,0.1)
          }
        }
        nimble_monitors[[paste0("alpha_", pred_p)]] <- paste0("alpha_", pred_p)

      }
    }else{
      if (!stringr::str_detect(pred_p, "s\\(") && pred_p != "1" && !stringr::str_detect(pred_p, "f\\(")) {
        if ((tolower(delay_link)=="log"|(tolower(delay_link)=="logit"))){
          nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:(D+1)) { \n")
          nimble_code <- paste0(nimble_code, aggregate_space, "  ", "alpha_", pred_p, "[d", aggregate_index, "] ~ ", nimble_priors$delay.re, "\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "} \n")
          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste0("alpha_", pred_p)]]<- matrix(stats::rnorm(nimble_constants$A*(nimble_constants$D+1),0,0.1), ncol=nimble_constants$A)
            }else{
              nimble_initial_values[[c]][[paste0("alpha_", pred_p)]]<- stats::rnorm((nimble_constants$D+1),0,0.1)
            }
          }
          nimble_monitors[[paste0("alpha_", pred_p)]] <- paste0("alpha_", pred_p)

        }else{
          nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:D) { \n")
          nimble_code <- paste0(nimble_code, aggregate_space, "  ", "alpha_", pred_p, "[d", aggregate_index, "] ~ ", nimble_priors$delay.re, "\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "} \n")
          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste0("alpha_", pred_p)]]<- matrix(stats::rnorm(nimble_constants$A*nimble_constants$D,0,0.1), ncol=nimble_constants$A)
            }else{
              nimble_initial_values[[c]][[paste0("alpha_", pred_p)]]<- stats::rnorm(nimble_constants$D,0,0.1)
            }
          }
          nimble_monitors[[paste0("alpha_", pred_p)]] <- paste0("alpha_", pred_p)

        }

      }
    }

  }

  # Add priors for spline coefficients

  if (!is.null(spline_terms_delay)) {
    # for(j in 1:length(spline_terms_delay)){
    #   nimble_code <- paste0(nimble_code, "  alpha_spline_", spline_terms_delay[j],"[1:N] <- g_spline_",spline_terms_delay[j],"[1:N,1:K_",spline_terms_delay[j],"]%*%kappa_alpha_", spline_terms_delay[j], "[1:K_",spline_terms_delay[j],"] \n")
    #   nimble_code <- paste0(nimble_code, "  for(k_", spline_terms_delay[j], " in 1:K_",spline_terms_delay[j],"){\n")
    #   nimble_code <- paste0(nimble_code,"    ", "kappa_alpha_",spline_terms_delay[j],"[k_",spline_terms_delay[j], "] ~ dnorm(0, 10)\n",sep='')
    #   nimble_code <- paste0(nimble_code, "  }\n")
    #
    # }
    for(j in 1:length(spline_terms_delay)){
      if(stringr::str_detect(tolower(delay_link), "survivor")){
        # nimble_constants[[paste("S_alpha_",pred_name,i,sep='')]]<-jagam_output$jags.data$S1[,(1+(i-1)*dim(jagam_output$jags.data$S1)[1]):(i*dim(jagam_output$jags.data$S1)[1])]
        nimble_code <- paste0(nimble_code, aggregate_space, "alpha_spline_",spline_terms_delay[j],"[1:N", aggregate_index,"] <- ","g_spline_",spline_terms_delay[j], "[1:N, 1:K_alpha_",pred_name,"]%*%kappa_alpha_", pred_name, "[1:K_alpha_",pred_name, aggregate_index, "] \n")
        nimble_code <- paste0(nimble_code, aggregate_space, "kappa_alpha_",spline_terms_delay[j],"[1:K_alpha_",spline_terms_delay[j] ,aggregate_index,"] ~ dmnorm(zeros[", "1:K_alpha_",spline_terms_delay[j],"],","omega_alpha_",spline_terms_delay[j],"[", "1:K_alpha_",spline_terms_delay[j], ", 1:K_alpha_",spline_terms_delay[j], aggregate_index,"])\n",sep='')
        nimble_constants$zeros<-rep(0,max(c(length(nimble_constants$zeros), nimble_constants[[paste("K_alpha_",spline_terms_delay[j],sep='')]])))
        # Penalty matrix for MVN distribution.
        nimble_code <- paste0(nimble_code, aggregate_space, "omega_alpha_",spline_terms_delay[j],"[1:K_alpha_",spline_terms_delay[j], ", 1:K_alpha_",spline_terms_delay[j], aggregate_index,"] <-",sep='')
        S_n<-(dim(jagam_output_delay[[j]]$jags.data$S1)[2]/dim(jagam_output_delay[[j]]$jags.data$S1)[1])
        for(i in 1:S_n){
          nimble_code <- paste0(nimble_code,"S_alpha_",spline_terms_delay[j],"[1:K_alpha_",spline_terms_delay[j], ", (1 + ",(i-1),"*K_alpha_",spline_terms_delay[j],"):(",i,"*K_alpha_",spline_terms_delay[j], ")]/sigma_alpha_", spline_terms_delay[j], i, aggregate_single, "^2 + ",sep='')
          if(i==S_n){
            nimble_code <- substr(nimble_code, 1, nchar(nimble_code) - 3)
            nimble_code <- paste0(nimble_code," \n",sep='')
          }
        }
        for(i in 1:S_n){
          nimble_code <- paste0(nimble_code, aggregate_space, "sigma_alpha_",spline_terms_delay[j],i, aggregate_single," ~ ", nimble_priors$delay.spline, " \n")
          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste("sigma_alpha_",spline_terms_delay[j],i,sep='')]]<- stats::runif(nimble_constants$A, 0, 1)

            }else{
              nimble_initial_values[[c]][[paste("sigma_alpha_",spline_terms_delay[j],i,sep='')]]<- stats::runif(1, 0, 1)
            }
          }
          nimble_monitors[[paste("sigma_alpha_",spline_terms_delay[j],i,sep='')]] <- paste("sigma_alpha_",spline_terms_delay[j],i,sep='')

        }
      }else{
        # nimble_constants[[paste("S_alpha_",pred_name,i,sep='')]]<-jagam_output$jags.data$S1[,(1+(i-1)*dim(jagam_output$jags.data$S1)[1]):(i*dim(jagam_output$jags.data$S1)[1])]
        nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:D){ \n")
        nimble_code <- paste0(nimble_code, aggregate_space,"  ", "alpha_spline_",spline_terms_delay[j],"[1:N, d", aggregate_index,"] <- ","g_spline_",spline_terms_delay[j], "[1:N, 1:K_alpha_",pred_name,"]%*%kappa_alpha_", pred_name, "[1:K_alpha_",pred_name,", d", aggregate_index, "] \n")
        nimble_code <- paste0(nimble_code, aggregate_space,"  ", "kappa_alpha_",spline_terms_delay[j],"[1:K_alpha_",spline_terms_delay[j], ", d" ,aggregate_index,"] ~ dmnorm(zeros[", "1:K_alpha_",spline_terms_delay[j],"],","omega_alpha_",spline_terms_delay[j],"[", "1:K_alpha_",spline_terms_delay[j], ", 1:K_alpha_",spline_terms_delay[j], ", d", aggregate_index,"])\n",sep='')
        nimble_constants$zeros<-rep(0,max(c(length(nimble_constants$zeros), nimble_constants[[paste("K_alpha_",spline_terms_delay[j],sep='')]])))
        # Penalty matrix for MVN distribution.
        nimble_code <- paste0(nimble_code, aggregate_space,"  ", "omega_alpha_",spline_terms_delay[j],"[1:K_alpha_",spline_terms_delay[j], ", 1:K_alpha_",spline_terms_delay[j],", d", aggregate_index,"] <-",sep='')
        S_n<-(dim(jagam_output_delay[[j]]$jags.data$S1)[2]/dim(jagam_output_delay[[j]]$jags.data$S1)[1])
        for(i in 1:S_n){
          nimble_code <- paste0(nimble_code,"S_alpha_",spline_terms_delay[j],"[1:K_alpha_",spline_terms_delay[j], ", (1 + ",(i-1),"*K_alpha_",spline_terms_delay[j],"):(",i,"*K_alpha_",spline_terms_delay[j], ")]/sigma_alpha_",spline_terms_delay[j],i,"[d",aggregate_index,"]^2 + ",sep='')
          if(i==S_n){
            nimble_code <- substr(nimble_code, 1, nchar(nimble_code) - 3)
            nimble_code <- paste0(nimble_code," \n",sep='')
          }
        }
        for(i in 1:S_n){
          nimble_code <- paste0(nimble_code, aggregate_space,"  ", "sigma_alpha_",spline_terms_delay[j],i,"[d",aggregate_index, "] ~ ", nimble_priors$delay.spline, " \n")
          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste("sigma_alpha_",spline_terms_delay[j],i,sep='')]]<- matrix(stats::runif(nimble_constants$A*nimble_constants$D, 0, 1),nrow=nimble_constants$D)

            }else{
              nimble_initial_values[[c]][[paste("sigma_alpha_",spline_terms_delay[j],i,sep='')]]<- stats::runif(nimble_constants$D, 0, 1)
            }
          }
          nimble_monitors[[paste("sigma_alpha_",spline_terms_delay[j],i,sep='')]] <- paste("sigma_alpha_",spline_terms_delay[j],i,sep='')

        }
        # End delay loop
        nimble_code <- paste0(nimble_code, aggregate_space," } \n")
      }
    }
  }

  # Add priors for random effects and structured terms
  if (length(random_effects_delay) > 0) {
    for (random_name in names(random_effects_delay)) {
      if(stringr::str_detect(tolower(delay_link), "survivor")){
        nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:",random_name,"_length) { \n")
        nimble_code <- paste0(nimble_code, aggregate_space,  "  ","v_", random_name, "[i", aggregate_index, "] ~ dnorm(0, tav_", random_name, aggregate_single ,")\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
        nimble_code <- paste0(nimble_code, aggregate_space,  "  ","tav_", random_name,  aggregate_single, " ~ ", nimble_priors$delay.iid, "\n")
      }else{
        if ((tolower(delay_link)=="log"|(tolower(delay_link)=="logit"))){
          nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:(D+1)) { \n")
        }else{
          nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:D) { \n")
        }
        nimble_code <- paste0(nimble_code, aggregate_space, "  ",  "for(i in 1:",random_name,"_length) { \n")
        nimble_code <- paste0(nimble_code, aggregate_space,  "  ", "  ","v_", random_name, "[i, d", aggregate_index, "] ~ dnorm(0, tav_", random_name, "[d", aggregate_index, "])\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
        nimble_code <- paste0(nimble_code, aggregate_space,  "  ","tav_", random_name, "[d", aggregate_index, "] ~ ", nimble_priors$delay.iid, "\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
      }

    }
  }

  if (length(structured_terms_delay) > 0) {
    for (term_name in names(structured_terms_delay)) {
      if(stringr::str_detect(tolower(delay_link), "survivor")){
        if (structured_terms_delay[[term_name]] == "rw1") {
          nimble_code <- paste0(nimble_code, aggregate_space, "v_rw1_", term_name, "[1", aggregate_index,"] ~ dnorm(0",", ", "tav_rw1_", term_name, aggregate_single,")\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "for(j_", term_name, " in 2:N){\n")
          nimble_code <- paste0(nimble_code, aggregate_space,"  ", "v_rw1_", term_name, "[j_", term_name, aggregate_index,"] ~ dnorm(v_rw1_", term_name,"[j_", term_name, "-1", aggregate_index, "]",", ", "tav_rw1_", term_name, aggregate_single, ")\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "tav_rw1_", term_name, aggregate_single, " ~ ", nimble_priors$delay.rw1, "\n")

          # } else if (structured_terms_delay[[term_name]] == "rw2") {
          #   nimble_code <- paste0(nimble_code, aggregate_space, "v_rw2_", term_name, "[1", aggregate_index,"] ~ dnorm(0, tav_rw2_", term_name, aggregate_single, ")\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "v_rw2_", term_name, "[2", aggregate_index,"] ~ dnorm(v_rw2_", term_name, "[1", aggregate_index,"], ", "tav_rw2_", term_name, aggregate_single, ")\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "for(j_", term_name, " in 3:N){\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "v_rw2_", term_name, "[j_", term_name, aggregate_index,"] ~ dnorm(2*v_rw2_", term_name,"[j_", term_name, "-1", aggregate_index,"] - v_rw2_", term_name, "[j_", term_name, "-2", aggregate_index,"]", ", ", "tav_rw2_", term_name, aggregate_single, ")\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "tav_rw2_", term_name, aggregate_single, " ~ ", nimble_priors$delay.rw2, "\n")
          #
          # } else if (structured_terms_delay[[term_name]] == "seasonal") {
          #   nimble_code <- paste0(nimble_code, "  v_seasonal_", term_name, "[1:N] ~ dPeriodic(tav_seasonal_", term_name, ")\n")
          #   nimble_code <- paste0(nimble_code, "  tav_seasonal_", term_name, " ~ T(dnorm(0,1),0,)\n")


          # }else if (structured_terms_delay[[term_name]] == "spatial") { #EDIT for spatial effect
          #   nimble_code <- paste0(nimble_code, aggregate_space, "v_spatial[1:N", aggregate_index,"] ~ dcar_normal(adj[], num[], tav_spatial",aggregate_single,")\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "tav_spatial",aggregate_single,," ~ T(dnorm(0,1),0,)\n")
          #
        }


      }else{
        if (structured_terms_delay[[term_name]] == "rw1") {
          if ((tolower(delay_link)=="log"|(tolower(delay_link)=="logit"))){
            nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:(D+1)) { \n")
          }else{
            nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:D) { \n")
          }
          nimble_code <- paste0(nimble_code, aggregate_space, "  ","v_rw1_", term_name, "[1, d", aggregate_index,"] ~ dnorm(0",", ", "tav_rw1_", term_name, "[d", aggregate_index, "])\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "  ","for(j_", term_name, " in 2:N){\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "  ","  ", "v_rw1_", term_name, "[j_", term_name, ", d", aggregate_index,"] ~ dnorm(v_rw1_", term_name,"[j_", term_name, "-1, d", aggregate_index, "]",", ", "tav_rw1_", term_name, "[d", aggregate_index, "])\n")
          nimble_code <- paste0(nimble_code, aggregate_space,"  ", "}\n")
          nimble_code <- paste0(nimble_code, aggregate_space,"  ", "tav_rw1_", term_name, "[d", aggregate_index, "] ~ ", nimble_priors$delay.rw1, "\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "}\n")

          # } else if (structured_terms_delay[[term_name]] == "rw2") {
          #   if ((tolower(delay_link)=="log"|(tolower(delay_link)=="logit"))){
          #     nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:(D+1)) { \n")
          #   }else{
          #     nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:D) { \n")
          #   }
          #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "v_rw2_", term_name, "[1, d", aggregate_index,"] ~ dnorm(0, tav_rw2_", term_name,"[d", aggregate_index, "])\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "v_rw2_", term_name, "[2, d", aggregate_index,"] ~ dnorm(v_rw2_", term_name, "[1, d", aggregate_index,"], ", "tav_rw2_", term_name,"[d", aggregate_index, "])\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "for(j_", term_name, " in 3:N){\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "  ", "v_rw2_", term_name, "[j_", term_name,", d", aggregate_index,"] ~ dnorm(2*v_rw2_", term_name,"[j_", term_name, "-1, d", aggregate_index,"] - v_rw2_", term_name, "[j_", term_name, "-2, d", aggregate_index,"]", ", ", "tav_rw2_", term_name,"[d", aggregate_index, "])\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "}\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "tav_rw2_", term_name,"[d", aggregate_index, "] ~ ", nimble_priors$delay.rw2, "\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
          #
          #   # } else if (structured_terms_delay[[term_name]] == "seasonal") {
          #   #   nimble_code <- paste0(nimble_code, "  v_seasonal_", term_name, "[1:N] ~ dPeriodic(tav_seasonal_", term_name, ")\n")
          #   #   nimble_code <- paste0(nimble_code, "  tav_seasonal_", term_name, " ~ T(dnorm(0,1),0,)\n")
          #
          # } else if (structured_terms_delay[[term_name]] == "spatial") { #EDIT for spatial effect
          #   nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:D) { \n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "v_spatial[1:N, d", aggregate_index,"] ~ dcar_normal(adj[], num[], tav_spatial[d", aggregate_index, "])\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "  ", "tav_spatial[d", aggregate_index, "] ~ T(dnorm(0,1),0,)\n")
          #   nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
          #
        }
      }

    }
  }

  if((tolower(family$delay)=="nb")|(tolower(family$delay)=="gdm")|(tolower(family$delay)=="dirichlet-multinomial")){
    # Add priors for family-specific parameters
    if((tolower(family$delay)=="gdm")){
      nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:D) {\n")
    }else{
      nimble_code <- paste0(nimble_code, aggregate_space, "for(d in 1:(D+1)) {\n")
    }
    nimble_code <- paste0(nimble_code, aggregate_space, "  ", "phi[d", aggregate_index,"] ~ ", nimble_priors$delay.dispersion, "\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
  }


  ### Priors for nested model

  if(nested){

    # EDIT need to add link choices for multiple viruses - modelling relative proportion!
    # Add prior distributions for the coefficients
    nimble_code <- paste0(nimble_code, aggregate_space, "# Priors nested\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "  gamma_0", aggregate_single , "~ ", nimble_priors$nested.intercept, "\n")
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][["gamma_0"]]<- stats::rnorm(nimble_constants$A,0,0.1)
      }else{
        nimble_initial_values[[c]][["gamma_0"]]<- stats::rnorm(1,0,0.1)
      }
    }
    nimble_monitors[["gamma_0"]] <- "gamma_0"



    # Add priors for regular predictors
    for (pred in predictors_nested) {
      if (!stringr::str_detect(pred, "s\\(") && pred != "1" && !stringr::str_detect(pred, "f\\(")) {
        nimble_code <- paste0(nimble_code, aggregate_space, "gamma_", pred, aggregate_single, " ~ ", nimble_priors$nested.re, "\n")
        # for(c in 1:nchains){
        #   if(aggregate){
        #     nimble_initial_values[[c]][[paste0("gamma_", pred)]]<- stats::rnorm(nimble_constants$A,0,0.1)
        #   }else{
        #     nimble_initial_values[[c]][[paste0("gamma_", pred)]]<- stats::rnorm(1,0,0.1)
        #   }
        # }
        # nimble_monitors[[paste0("gamma_", pred)]] <- paste0("gamma_", pred)

      }
    }

    # Add priors for spline coefficients
    if (!is.null(spline_terms_nested)) {
      nimble_code <- paste0(nimble_code, aggregate_space, "# Priors for totals spline \n")
      for(j in 1:length(spline_terms_nested)){
        # nimble_constants[[paste("S_gamma_",pred_name,i,sep='')]]<-jagam_output$jags.data$S1[,(1+(i-1)*dim(jagam_output$jags.data$S1)[1]):(i*dim(jagam_output$jags.data$S1)[1])]
        nimble_code <- paste0(nimble_code, aggregate_space, "gamma_spline_",spline_terms_nested[j],"[1:N", aggregate_index, "] <- ","f_spline_",spline_terms_nested[j], "[1:N, 1:K_gamma_",pred_name,"]%*%kappa_gamma_", pred_name, "[1:K_gamma_",pred_name, aggregate_index, "] \n")
        nimble_code <- paste0(nimble_code, aggregate_space, "kappa_gamma_",spline_terms_nested[j],"[", "1:K_gamma_",spline_terms_nested[j], aggregate_index, "] ~ dmnorm(zeros[", "1:K_gamma_",spline_terms_nested[j],"],","omega_gamma_",spline_terms_nested[j],"[", "1:K_gamma_",spline_terms_nested[j], ", 1:K_gamma_",spline_terms_nested[j], aggregate_index, "])\n",sep='')
        nimble_constants$zeros<-rep(0,max(c(length(nimble_constants$zeros), nimble_constants[[paste("K_gamma_",spline_terms_nested[j],sep='')]])))
        # Penalty matrix for MVN distribution.
        nimble_code <- paste0(nimble_code, aggregate_space, "omega_gamma_",spline_terms_nested[j],"[", "1:K_gamma_",spline_terms_nested[j], ", 1:K_gamma_",spline_terms_nested[j], aggregate_index, "] <-",sep='')
        S_n<-(dim(jagam_output_nested[[j]]$jags.data$S1)[2]/dim(jagam_output_nested[[j]]$jags.data$S1)[1])
        for(i in 1:S_n){
          nimble_code <- paste0(nimble_code,"S_gamma_",spline_terms_nested[j],"[", "1:K_gamma_",spline_terms_nested[j], ", (1 + ",(i-1),"*K_gamma_",spline_terms_nested[j],"):(",i,"*K_gamma_",spline_terms_nested[j], ")]/sigma_gamma_",spline_terms_nested[j],i,aggregate_single,"^2 + ",sep='')
          if(i==S_n){
            nimble_code <- substr(nimble_code, 1, nchar(nimble_code) - 3)
            nimble_code <- paste0(nimble_code," \n",sep='')
          }
        }
        for(i in 1:S_n){
          nimble_code <- paste0(nimble_code, aggregate_space, "sigma_gamma_", spline_terms_nested[j],i,aggregate_single," ~ ", nimble_priors$nested.spline, " \n")
          for(c in 1:nchains){
            if(aggregate){
              nimble_initial_values[[c]][[paste("sigma_gamma_", spline_terms_nested[j],i,sep='')]]<- stats::runif(nimble_constants$A, 0, 1)
            }else{
              nimble_initial_values[[c]][[paste("sigma_gamma_", spline_terms_nested[j],i,sep='')]]<- stats::runif(1, 0, 1)
            }
            nimble_monitors[[paste("sigma_gamma_", spline_terms_nested[j],i,sep='')]] <- paste("sigma_gamma_", spline_terms_nested[j],i,sep='')
          }
        }
      }
    }

    # Add priors for random effects and structured terms
    if (length(random_effects_nested) > 0) {
      for (random_name in names(random_effects_nested)) {
        nimble_code <- paste0(nimble_code, aggregate_space, "for(i in 1:",random_name,"_length){\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "  ", "w_", random_name, "[i",aggregate_index,"] ~ dnorm(0, taw_", random_name, aggregate_single, ")\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
        nimble_code <- paste0(nimble_code, aggregate_space, "taw_", random_name, aggregate_single, " ~ ", nimble_priors$nested.iid, "\n")
      }
    }

    if (length(structured_terms_nested) > 0) {
      for (term_name in names(structured_terms_nested)) {
        if (structured_terms_nested[[term_name]] == "rw1") {
          nimble_code <- paste0(nimble_code, aggregate_space, "w_rw1_", term_name, "[1", aggregate_index, "] ~ dnorm(0",", ", "taw_rw1_", term_name, ")\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "for(j_", term_name, " in 2:N){\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "  ", "w_rw1_", term_name, "[j_", term_name, aggregate_index, "] ~ dnorm(w_rw1_", term_name,"[j_", term_name, "-1]",", ", "taw_rw1_", term_name, ")\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "taw_rw1_", term_name, aggregate_single, " ~ ", nimble_priors$nested.rw1, "\n")
        } else if (structured_terms_nested[[term_name]] == "rw2") {
          nimble_code <- paste0(nimble_code, aggregate_space, "w_rw2_", term_name, "[1", aggregate_index,"] ~ dnorm(0, taw_rw2_", term_name, aggregate_single, ")\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "w_rw2_", term_name, "[2", aggregate_index,"] ~ dnorm(w_rw2_", term_name, "[1",aggregate_index,"], ", "taw_rw2_", term_name, aggregate_single, ")\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "for(j_", term_name, " in 3:N){\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "    ", "w_rw2_", term_name, "[j_", term_name, aggregate_index , "] ~ dnorm(2*w_rw2_", term_name,"[j_", term_name, "-1", aggregate_index,"] - w_rw2_", term_name, "[j_", term_name, "-2", aggregate_index,"]", ", ", "taw_rw2_", term_name, aggregate_single, ")\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "}\n")
          nimble_code <- paste0(nimble_code, aggregate_space, "taw_rw2_", term_name, aggregate_single, " ~ ", nimble_priors$nested.rw2, "\n")
          # } else if (structured_terms_nested[[term_name]] == "seasonal") {
          #   nimble_code <- paste0(nimble_code, "  w_seasonal_", term_name, "[1:N] ~ dPeriodic(taw_seasonal_", term_name, ")\n")
          #   nimble_code <- paste0(nimble_code, "  taw_seasonal_", term_name, " ~ T(dnorm(0,1),0,)\n")
        }
        # else if (structured_terms_nested[[term_name]] == "spatial") {
        #   nimble_code <- paste0(nimble_code, aggregate_space, "w_spatial[1:N", aggregate_index,"] ~ dcar_normal(adj[], num[], taw_spatial",aggregate_single,")\n")
        #   nimble_code <- paste0(nimble_code, aggregate_space, "taw_spatial",aggregate_single," ~ T(dnorm(0,1),0,)\n")
        # }
      }
    }
    # Initial values for delta and omega
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][["delta"]]<- stats::rnorm(nimble_constants$A,0,1)
        nimble_initial_values[[c]][["omega"]]<- stats::runif(nimble_constants$A,0,0.1)

      }else{
        nimble_initial_values[[c]][["delta"]]<- stats::rnorm(1,0,1)
        nimble_initial_values[[c]][["omega"]]<- stats::runif(1,0,0.1)
      }
    }


    # Prior for scale of shared parameter
    nimble_code <- paste0(nimble_code, aggregate_space, "delta", aggregate_single, " ~ ", nimble_priors$nested.scale, "\n")

    # Prior for censoring slope
    nimble_code <- paste0(nimble_code, aggregate_space, "omega", aggregate_single, " ~ ", nimble_priors$censored.slope, "\n")

    # Add priors for family-specific parameters
    nimble_code <- paste0(nimble_code, aggregate_space, "upsilon", aggregate_single, " ~ ", nimble_priors$nested.dispersion, "\n")
    nimble_code <- paste0(nimble_code, aggregate_space, "chi", aggregate_single, " ~ ", nimble_priors$censored.dispersion, "\n")

  }



  # End of nimble code:

  # If aggregate is TRUE close additional for loop
  if(aggregate){
    nimble_code<-paste0(nimble_code, "  }\n", sep='')
  }

  # Close "nimble_code({" function
  nimble_code <- paste0(nimble_code, "})")

  ############## Lists for nimble:


  # If totals response variable is not in the data calculate from partial response variable data:
  if(!conditional){
    if(is.null(data[[response]])){
      if(aggregate){
        data[[response]]<-apply(data[[response_partial]],c(1,3),sum)

      }else{
        data[[response]]<-apply(data[[response_partial]],1,sum)
      }
    }

  }


  if(is.null(model_window)){
    if(!conditional){
      if(aggregate){
        nimble_data[[response]]<-data[[response]][1:nimble_constants$W,]

      }else{
        nimble_data[[response]]<-data[[response]][1:nimble_constants$W]

      }
    }

    if(aggregate){
      nimble_data[[response_partial]]<-data[[response_partial]][1:nimble_constants$W,1:D_index,]

    }else{
      nimble_data[[response_partial]]<-data[[response_partial]][1:nimble_constants$W,1:D_index]

    }
  }else{
    if(aggregate){
      if(!conditional){
        nimble_data[[response]]<-data[[response]][(1+dim(data[[response]])[1]-nimble_constants$W-forecast):(dim(data[[response]])[1]-forecast),]
      }
      nimble_data[[response_partial]]<-data[[response_partial]][(1+dim(data[[response]])[1]-nimble_constants$W-forecast):(dim(data[[response]])[1]-forecast),1:D_index,]


    }else{
      if(!conditional){
        nimble_data[[response]]<-data[[response]][(1+length(data[[response]])-nimble_constants$W-forecast):(length(data[[response]])-forecast)]
      }
      nimble_data[[response_partial]]<-data[[response_partial]][(1+length(data[[response]])-nimble_constants$W-forecast):(length(data[[response]])-forecast),1:D_index]
    }
  }



  if(nested){
    if(aggregate){
      nimble_data[[response_nested]]<-data[[response_nested]][1:nimble_constants$W,]

    }else{
      nimble_data[[response_nested]]<-data[[response_nested]][1:nimble_constants$W]

    }
    observed_nested<-nimble_data[[response_nested]]
    if(aggregate){
      observed_nested[(nimble_constants$C+1):(nimble_constants$W),]<-NA

    }else{
      observed_nested[(nimble_constants$C+1):(nimble_constants$W)]<-NA

    }
    nimble_data[[paste0(response_nested,"_corrected")]]<-observed_nested
  }

  # Missing data initial values for nimble
  #EDIT initial values for forecasting period? two time loops?
  #EDIT if for forecast not with NA in data?

  if(aggregate){
    # Don't need to set initial values if partial counts aren't being sampled
    if(!(tolower(family$delay)=="gdm"|tolower(family$delay)=="nb"|tolower(family$delay)=="poisson")){
      if(any(is.na(nimble_data[[response_partial]]))){
        for(c in 1:nchains){
          nimble_initial_values[[c]][[response_partial]] <- array(NA, dim=dim(nimble_data[[response_partial]]))
          if(!conditional){
            nimble_initial_values[[c]][[response]] <- array(NA, dim=dim(nimble_data[[response]]))
          }
          for(a in 1:nimble_constants$A){
            for(t in 1:nimble_constants$W){
              for(d in 1:(nimble_constants$D+1)){
                # Initial values for unobserved response_delay (delay).
                if(is.na(nimble_data[[response_partial]][t,d,a])){
                  nimble_initial_values[[c]][[response_partial]][t,d,a]<- stats::rpois(1, stats::median(nimble_data[[response_partial]][,d,a],na.rm=T))
                }
              }
              # Initial values for unobserved response (totals).
              if(!conditional){
                if(is.na(nimble_data[[response]][t,a])){
                  nimble_initial_values[[c]][[response]][t,a] <- sum(nimble_initial_values[[c]][[response_partial]][t,,a],na.rm = TRUE) +
                    sum(nimble_data[[response_partial]][t,,a], na.rm = TRUE)
                }
              }
            }
          }
        }
      }
    }else if(!conditional){
      # Initial values for unobserved response (totals) for gdm/nb/poisson model which isn't a conditional indepedence model.
      if(any(is.na(nimble_data[[response]]))){
        for(c in 1:nchains){
          nimble_initial_values[[c]][[response]] <- array(NA, dim=dim(nimble_data[[response]]))
          for(s in 1:nimble_constants$A){
            nimble_initial_values[[c]][[response]][which(is.na(nimble_data[[response]][,s])),s]<-#apply(nimble_initial_values[[c]][[response_partial]][which(is.na(nimble_data[[response]][,s])),,s],1,sum,na.rm = TRUE) +
              apply(nimble_data[[response_partial]][which(is.na(nimble_data[[response]][,s])),,s],1,sum,na.rm = TRUE) +
              stats::rpois(length(which(is.na(nimble_data[[response]][,s]))),10) +
              stats::rpois(length(which(is.na(nimble_data[[response]][,s]))),stats::median(nimble_data[[response]][,s]-rowSums(nimble_data[[response_partial]][,,s]),na.rm=T))
          }
        }
      }
    }

  }else{
    # Don't need to set initial values if partial counts aren't being sampled.
    if(!(tolower(family$delay)=="gdm"|tolower(family$delay)=="nb"|tolower(family$delay)=="poisson")){
      if(any(is.na(nimble_data[[response_partial]]))){
        for(c in 1:nchains){

          nimble_initial_values[[c]][[response_partial]] <- array(NA, dim=dim(nimble_data[[response_partial]]))
          if(!conditional){
            nimble_initial_values[[c]][[response]] <- rep(NA, length(nimble_data[[response]]))
          }
          for(t in 1:nimble_constants$N){
            for(d in 1:(nimble_constants$D+1)){
              # Initial values for unobserved response_delay (delay).
              if(is.na(nimble_data[[response_partial]][t,d])){
                nimble_initial_values[[c]][[response_partial]][t,d]<- stats::rpois(1, stats::median(nimble_data[[response_partial]][,d],na.rm=T))
              }
            }
            # Initial values for unobserved response (totals) for gdm/nb/poisson model which isn't a conditional indepedence model.
            if(!conditional){
              if(is.na(nimble_data[[response]][t])){
                nimble_initial_values[[c]][[response]][t] <- sum(nimble_initial_values[[c]][[response_partial]][t,], na.rm = TRUE) +
                  sum(nimble_data[[response_partial]][t,], na.rm = TRUE)
              }
            }
          }
        }
      }
    }else if(!conditional){
      if(any(is.na(nimble_data[[response]]))){
        for(c in 1:nchains){

          nimble_initial_values[[c]][[response]] <- rep(NA, length(nimble_data[[response]]))
          nimble_initial_values[[c]][[response]][which(is.na(nimble_data[[response]]))]<-#apply(nimble_initial_values[[c]][[response_partial]][which(is.na(nimble_data[[response]])),],1,sum,na.rm = TRUE)
            apply(nimble_data[[response_partial]][which(is.na(nimble_data[[response]])),],1,sum,na.rm = TRUE) +
            stats::rpois(length(which(is.na(nimble_data[[response]]))),10)+
            stats::rpois(length(which(is.na(nimble_data[[response]]))),stats::median(nimble_data[[response]]-rowSums(nimble_data[[response_partial]]),na.rm=T))

        }
      }
    }
  }

  # Add initial values for corrected nested response
  if(nested){
    for(c in 1:nchains){
      if(aggregate){
        nimble_initial_values[[c]][[paste0(response_nested,"_corrected")]] <- array(NA, dim=dim(nimble_data[[paste0(response_nested,"_corrected")]]))

        for(a in 1:nimble_constants$A){
          for(i in (nimble_constants$C+1):nimble_constants$W){
            nimble_initial_values[[c]][[paste0(response_nested,"_corrected")]][i,a] <-
              ceiling(mean(c(nimble_data[[response_nested]][i,a],sum(nimble_data[[response]][i,a],nimble_initial_values[[c]][[response]][i,a],na.rm=T))))
          }
        }
        # if(nimble_constants$N>nimble_constants$W){
        #   for(a in 1:nimble_constants$A){
        #     for(i in (nimble_constants$W+1):nimble_constants$N){
        #       nimble_initial_values[[c]][[paste0(response_nested,"_corrected")]][i,a] <- stats::runif(1, min=0, max=nimble_initial_values[[c]][[response]][i,a])
        #     }
        #   }
        # }
      }else{
        nimble_initial_values[[c]][[paste0(response_nested,"_corrected")]] <- rep(NA, length(nimble_data[[paste0(response_nested,"_corrected")]]))

        for(i in (nimble_constants$C+1):nimble_constants$W){
          nimble_initial_values[[c]][[paste0(response_nested,"_corrected")]][i] <-
            ceiling(mean(c(nimble_data[[response_nested]][i],sum(nimble_data[[response]][i],nimble_initial_values[[c]][[response]][i],na.rm=T))))
        }

        # if(nimble_constants$N>nimble_constants$W){
        #   for(i in (nimble_constants$W+1):nimble_constants$N){
        #     nimble_initial_values[[c]][[paste0(response_nested,"_corrected")]][i] <- stats::runif(1, min=0, max=nimble_initial_values[[c]][[response]][i])
        #   }
        #
        # }
      }
    }
  }


  # Constants for nimble
  if(tolower(family$delay)=="gdm"|tolower(family$delay)=="nb"|tolower(family$delay)=="poisson"){
    if(aggregate){
      obs_index<-matrix(rep(nimble_constants$W, D_index*nimble_constants$A), ncol=nimble_constants$A)  # Nowcasting and forecasting time steps
      for(d in 1:D_index){
        for(s in 1:nimble_constants$A){
          if(any(is.na(data[[response_partial]][,d,s]))){
            obs_index[d,s]<-which(is.na(data[[response_partial]][,d,s])==TRUE)[1]-1
          }
        }
      }
    }else{
      obs_index<-rep(nimble_constants$W, D_index)  # Nowcasting and forecasting time steps
      for(d in 1:D_index){
        if(any(is.na(data[[response_partial]][,d]))){
          obs_index[d]<-which(is.na(data[[response_partial]][,d])==TRUE)[1]-1
        }
      }
    }
    nimble_constants$obs_index<-obs_index

  }


  # Constant for model window length not needed if delay model
  if(conditional){
    nimble_constants$W<-NULL
  }



  return(list(nimble_code = nimble_code, nimble_initial_values = nimble_initial_values, #nimble_dimensions = nimble_dimensions,
              nimble_constants = nimble_constants, nimble_data=nimble_data, nimble_monitors=nimble_monitors))
}

