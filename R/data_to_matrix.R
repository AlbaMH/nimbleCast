#' Convert individual patient data to a matrix.
#'
#' @param data.individual Data set formatted with individual patient level data.
#' @param window Number of weeks from the most recent date you want to keep.
#' @param D_max Number of delays in weeks you want to be returned. Output will have D_max+1 delay dimensions where D_max+1 is all remainder counts after D_max.
#' @param nested_col Name of the column related to a potential nested structure you wish to model. Can be set to "default.covid" if COVID-19 test results are in columns "AN_SARS2" & "PCR_SARS2".
#' @param bins_age A numeric vector of bins you wish to divide age into e.g. c(0,25,50,75,100,125), or can be set to "SI-PNI", "5 years" or "10 year".
#' @param date_onset The column name corresponding to the date of onset.
#' @param date_report The column name corresponding to the date of reporting.
#' @param age_col The column name corresponding to the date of birth.
#' @param region_col The column name corresponding to the spatial region.
#' @param date_format The format of the dates in data.individual, default is "%d/%m/%Y".
#' @param forecast Integer of weeks of forecasting you want to include (to be passed to nimbleCast).
#' @param use.epiweek If TRUE, it uses the CDC epiweek definition where the week starts on Sunday, if FALSE it the week ends at the weekday of the last record date.
#' @param trim.data Integer such that nowcasting is performed up to trim.data weeks prior to the maximum date in the data.
#' @param testing Default is FALSE which implys data is most recent data available and nowcast date is most recent date_report. Set to TRUE if you wish to nowcast up to the most recent date_onset. If set as date in format yyyy-mm-dd that will be the nowcast date. If not FALSE arrays with all available data will also be returned.
#'
#' @import dplyr
#'
#' @return List of data arrays.
#' @export
#'
#' @examples
#' Sari to matrix
data_to_matrix <- function(data.individual,
                           window = NULL,
                           D_max= NULL,  # output matrix will have D_max+1 delays where D_max+1 is all cases reported after delay D_max
                           nested_col = c(FALSE, "default.covid", nested_col),
                           bins_age = c(NULL,"SI-PNI", "10 years", "5 years", bins_age),
                           date_onset,
                           date_report,
                           age_col=NULL, # Allow for age as well as DofB?
                           region_col=NULL,
                           date_format=c("%d/%m/%Y"),
                           forecast=0,
                           use.epiweek = FALSE,
                           trim.data,
                           testing=FALSE){

  # Clean data set dates
  # if nested is not false create nesed matrix too
  if(tolower(nested_col)=="default.covid"){
    col_ant<-which(colnames(data.individual)=="AN_SARS2")
    col_pcr<-which(colnames(data.individual)=="PCR_SARS2")
    which_nested<-c(col_ant,col_pcr)
  }else if(!(nested_col==FALSE)){
    which_nested<-which(colnames(data.individual)==nested_col)
  }else{
    which_nested<-NULL
  }

  dataset <- data.individual[,c(which_nested, which(colnames(data.individual)==date_onset),which(colnames(data.individual)==region_col),which(colnames(data.individual)==date_report), which(colnames(data.individual)==age_col))]
  dataset[[date_onset]]<-as.Date((dataset[[date_onset]]),format=date_format)
  dataset[[date_report]]<-as.Date((dataset[[date_report]]),format=date_format)
  if(!is.null(age_col)){
    dataset[[age_col]]<-as.Date((dataset[[age_col]]),format=date_format)

  }

  if(tolower(nested_col)=="default.covid"){
    dataset[[nested_col]]<-as.numeric((dataset$PCR_SARS2==1)|(dataset$AN_SARS2==1))
    dataset[is.na(dataset[[nested_col]]),nested_col]<-0
    dataset$PCR_SARS2<-NULL
    dataset$AN_SARS2<-NULL
  }

  # Set bins for ages:
  if(!is.null(age_col)){
    if(!is.numeric(bins_age)){
      if(missing(bins_age) | is.null(bins_age)){
        bins_age<-c(0,5,12,18,30,seq(40,90,by=10),130)
        labels_age<-1:(length(bins_age)-1)
        print(paste0("Bins age as in SI-PNI: ",
                     paste0(bins_age[bins_age < 90], collapse = " "), " ",  "90+"))
      } else if(bins_age == "SI-PNI"){
        bins_age<-c(0,5,12,18,30,seq(40,90,by=10),130)
        labels_age<-1:(length(bins_age)-1)
        print(paste0("Bins age as in SI-PNI: ",
                     paste0(bins_age[bins_age < 90], collapse = " "), " ",  " 90+"
                ))
      } else if(bins_age == "10 years"){
        bins_age<-c(seq(0,90, by = 10),130)
        labels_age<-1:(length(bins_age)-1)
        print(paste0("Bins age in 10 years: ",
                     paste0(bins_age[bins_age < 90], collapse = " "), " ",  "90+"))
      } else if(bins_age == "5 years"){
        bins_age<-c(seq(0,90, by = 5),130)
        labels_age<-1:(length(bins_age)-1)
        print(paste0("Bins age in 5 years: ",
                     paste0(bins_age[bins_age < 90], collapse = " "), " ", "90+"))
      }
      else {
        stop("Age bins are not options of 'SI-PNI', '10 years', '5 years' or numeric vector!")
      }
    } else {
      bins_age<-bins_age
      labels_age<-1:(length(bins_age)-1)
      print(paste0("Using bins ages given: ",
              paste0(bins_age[bins_age < bins_age[length(bins_age) - 1]], collapse =" "), " ",
              bins_age[length(bins_age) - 1], "+"))
    }
  }else{
    bins_age<-c(0,130)
    labels_age<-1
    print("No age_col given so data will be returned withoout binning by age groups.")
  }


  ## Last digitalisation date considered
  if(missing(trim.data)){
    trim.data <-  0
    print("Using default, no trimming out of the data")
  } else {
    print(paste0("Trimming out ",
                 trim.data ,
                 " weeks of data"),
          call. = T)
  }

  ## Transforming trim.data into days
  trim.data.w <- 7*trim.data

  if(testing==FALSE){

    ## Maximum date to be considered on the estimation
    DT_max <- max(dataset |>
                    dplyr::pull(var = {{date_report}}),
                  na.rm = T) - trim.data.w
    if(!is.null(D_max)){
      if(DT_max>(max(dataset |>
                     dplyr::pull(var = {{date_onset}}),
                     na.rm = T)+D_max*7+1) ){
        DT_max<-max(dataset |>
                      dplyr::pull(var = {{date_onset}}),
                    na.rm = T)- trim.data.w + D_max*7 + 1
        print(paste0("Maximum date_report is over D_max weeks greater than maximum date_onset. Setting nowcast date as ",DT_max,"."))

      }else{
        print(paste0("Setting nowcast date as ",DT_max," (maximum date_report in data)."))

      }
    }else{
      if(DT_max>(max(dataset |>
                     dplyr::pull(var = {{date_onset}}),
                     na.rm = T)) ){
        print(paste0("Setting nowcast date as ",DT_max," (maximum date_report in data). Please note this is ",ceiling(as.numeric((DT_max-(max(dataset |>
                                                                                                                                                dplyr::pull(var = {{date_onset}}),
                                                                                                                                              na.rm = T))))/7) ," weeks greater than the most recent onset date."))
      }else{
        print(paste0("Setting nowcast date as ",DT_max," (maximum date_report in data)."))

      }
    }

  }else{

    ## Maximum date to be considered on the estimation
    if(testing==TRUE){
      ## If testing=TRUE data is historic and date reported may be greater than date onset:
      DT_max <- max(dataset |>
                      dplyr::pull(var = {{date_onset}}),
                    na.rm = T)
      print(paste0("Setting nowcast date as ",DT_max," (maximum date_onset in data as testing=TRUE)."))

    }else{
      DT_max <- as.Date(testing)
      print(paste0("Setting nowcast date as ",DT_max," as set by testing."))

    }

  }

  ## Transforming window into days
  if(is.null(window)){
    window.data.w<- max(dataset |>
                          dplyr::pull(var = {{date_report}}),
                        na.rm = T) - min(dataset |>
                                           dplyr::pull(var = {{date_onset}}),
                                         na.rm = T)
  }else{
    window.data.w <- 7*window

  }

  ## Last day of the week for the digitisation date calculation
  DT_max_diadasemana <- as.integer(format(DT_max, "%w"))

  #  Notice that if max recording date is Saturday (DT_max_diadasemana = 6) then the week is complete,
  # and the epiweek is on course. Otherwise some data must be ignored

  ## Ignore data after the last Sunday of recording (Sunday as too)
  aux.trimming.date = ifelse( use.epiweek | DT_max_diadasemana == 6, DT_max_diadasemana + 1, 0)


  if(!(testing==FALSE)){
    ## Calculate historic totals if testing isn't FALSE
    data_true <- dataset |>
      dplyr::rename(date_report = {{date_report}},
                    date_onset = {{date_onset}},
                    age_col = {{age_col}},
                    nested_col ={{nested_col}},
                    region_col ={{region_col}}) |>
      dplyr::filter(date_onset <= DT_max - aux.trimming.date + forecast*7,
                    date_onset>=(DT_max - aux.trimming.date - window.data.w + 1)) |>
      tidyr::drop_na(age_col) |>
      dplyr::mutate(
        ## Onset date
        # Moving the date to sunday
        DT.sun.aux = as.integer(format(date_onset, "%w")),
        ## Altering the date for the first day of the week
        dt.aux = date_onset -
          # Last recording date (DT_max_diadasemana) is the last day of the new week format
          DT.sun.aux +
          ifelse( use.epiweek, 0, DT_max_diadasemana+1 -
                    ifelse(DT_max_diadasemana+1>DT.sun.aux,7, 0)
          ),
        date_onset = dt.aux - ifelse( date_onset < dt.aux, 7, 0),
        # Recording date
        DT.sun.aux.rep = as.integer(format(date_report, "%w")),
        ## Altering the date for the first day of the week
        dt.aux.rep = date_report -
          # Last recording date (DT_max_diadasemana) is the last day of the new week format
          DT.sun.aux.rep +
          ifelse( use.epiweek, 0, DT_max_diadasemana+1 -
                    ifelse(DT_max_diadasemana+1 > DT.sun.aux.rep,7, 0)
          ),
        date_report = dt.aux.rep - ifelse( date_report < dt.aux.rep, 7, 0),
        Delay = as.numeric(date_report - date_onset) / 7,

      ) |>
      dplyr::select(-dt.aux,-dt.aux.rep, -DT.sun.aux,-DT.sun.aux.rep) |>
      dplyr::filter(Delay >= 0)


    ## Minimum date of onset
    first_date_true <- min(data_true |>
                             dplyr::pull(var = date_onset))

    ## Calculate onset weeks
    data_true<-dplyr::mutate(data_true, onset_week = (as.numeric((date_onset-first_date_true)/7)+1))

    # Create age col and remove any data greater than maximum age bin
    if(!is.null(age_col)){
      data_true <- data_true |>
        dplyr::mutate(age = floor(as.numeric(date_onset-age_col)/365.24),
                      fx_etaria = cut(age,  breaks = bins_age, labels = labels_age, right = F))|>
        dplyr::filter(age <= max(bins_age)) |>
        tidyr::drop_na(fx_etaria)
      # Calculate number of cases for each week
      # Maximum regions
      if(!is.null(region_col)){
        sum_cases_true<-data_true|>dplyr::group_by(date_report,date_onset,onset_week,Delay,fx_etaria,region_col)|>dplyr::summarise(nested_col=sum(nested_col),cases=dplyr::n())
        sum_cases_true<-sum_cases_true|> dplyr::mutate(region_col=as.factor(region_col))|> dplyr::mutate(region=as.numeric(region_col))
        region_length_true<-max((sum_cases_true$region))
        region_names<-levels(sum_cases_true$region_col)
      }else{
        sum_cases_true<-data_true|>dplyr::group_by(date_report,date_onset,onset_week,Delay,fx_etaria)|>dplyr::summarise(nested_col=sum(nested_col),cases=dplyr::n())
        region_length_true=1
      }
      Age_length_true<-length(levels(sum_cases_true$fx_etaria))

    }else{
      # Calculate number of cases for each week
      if(!is.null(region_col)){
        sum_cases_true<-data_true|>dplyr::group_by(date_report,date_onset,onset_week,Delay,region_col)|>dplyr::summarise(nested_col=sum(nested_col),cases=dplyr::n())
        sum_cases_true<-sum_cases_true|> dplyr::mutate(region_col=as.factor(region_col))|> dplyr::mutate(region=as.numeric(region_col))
        region_length_true<-max((sum_cases_true$region))
        region_names<-levels(sum_cases_true$region_col)
      }else{
        sum_cases_true<-data_true|>dplyr::group_by(date_report,date_onset,onset_week,Delay)|>dplyr::summarise(nested_col=sum(nested_col),cases=dplyr::n())
        region_length_true=1
        region_names<-NULL
      }
      Age_length_true<-1
      sum_cases_true<-sum_cases_true|>dplyr::mutate(fx_etaria=Age_length_true)
    }
    if(!is.null(region_col)){
      sum_cases_true<-sum_cases_true[order(sum_cases_true$region_col),]

    }
    sum_cases_true<-sum_cases_true[order(sum_cases_true$onset_week),]

    # Maximum onset week
    N_raw_true <- max(sum_cases_true$onset_week)

    # Maximum delay in data
    D_raw_true<-max(as.numeric(sum_cases_true$Delay),0)
    totals_list_true<-list()
    nested_list_true<-list()

    for(a in 1:Age_length_true){
      # Number of delays in data
      totals_list_unordered_true<-tidyr::pivot_wider(dplyr::select(dplyr::filter(sum_cases_true, fx_etaria==a),c(fx_etaria,date_report,date_onset,Delay,cases,onset_week,region_col)) ,id_cols=c(onset_week,region_col),names_from = Delay,values_from = cases)
      totals_list_unordered_true[is.na(totals_list_unordered_true)]<-0
      if(!is.null(region_col)){
        totals_list_unordered_true<-totals_list_unordered_true[order(totals_list_unordered_true$region_col),]

      }
      totals_list_unordered_true<-totals_list_unordered_true[order(totals_list_unordered_true$onset_week),]

      # Order delays
      order_delay_true<-totals_list_unordered_true[,-c(1:length(c(date_onset,region_col)))]
      order_delay_true<-order_delay_true[,order(as.integer(colnames(order_delay_true)))]
      totals_list_true_ordered1<-cbind(totals_list_unordered_true[,c(1:length(c(date_onset,region_col)))],order_delay_true)

      # Order attributes in rows:
      totals_list_true_ordered2<-totals_list_true_ordered1[,1:dim(totals_list_true_ordered1)[2]]|> dplyr::mutate(id=paste(onset_week,"-", region_col, sep=''))

      # Create blank data frame for all regions and weeks.
      if(is.null(region_col)){
        totals_list_true[[a]] <- cbind(tibble::tibble(onset_week=sort(rep((1:N_raw_true)))),
                                       matrix(0, ncol=D_raw_true+1,nrow=N_raw_true))|>dplyr::mutate(check = NA)|>dplyr::mutate(id=paste(onset_week,sep=''))

      }else{
        totals_list_true[[a]] <- cbind(tibble::tibble(onset_week=sort(rep((1:N_raw_true),region_length_true)), region_col=rep(region_names,N_raw_true)),
                                       matrix(0, ncol=D_raw_true+1,nrow=N_raw_true*region_length_true))|>dplyr::mutate(check = NA)|>dplyr::mutate(id=paste(onset_week,"-",region_col,sep=''))

      }

      # Fill in available data:
      totals_list_true[[a]][totals_list_true[[a]]$id %in% totals_list_true_ordered2$id, c(rep(FALSE,length(c(date_onset, region_col))),0:D_raw_true %in% as.numeric(colnames(order_delay_true)),TRUE,FALSE)] <- totals_list_true_ordered2[,(length(c(date_onset,region_col))+1):(dim(totals_list_true_ordered2)[2])]

      # Create D_max delays
      if(!is.null(D_max)){
        totals_list_true[[a]]<-cbind(totals_list_true[[a]][,(length(c(date_onset, region_col))+1):(D_max+length(c(date_onset, region_col)))],apply(totals_list_true[[a]][,(D_max+length(c(date_onset, region_col))+1):(dim(totals_list_true[[a]])[2]-2)],1,sum, na.rm=T))
        colnames(totals_list_true[[a]])[dim(totals_list_true[[a]])[2]]<-paste(D_max+1)
        D_dim_true<-D_max
      }else{
        D_dim_true<-D_raw_true
      }

      # Make a dimension for region if including region: array[REGION, TIME, DELAY]
      if(!is.null(region_col)){
        totals_list_true[[a]]<-totals_list_true[[a]][,1:(D_dim_true+1)]|>as.matrix()|>array(dim=c(region_length_true,N_raw_true,D_dim_true+1))
      }

      if(!is.null(nested_col)){
        if(is.null(region_col)){
          nested_unordered_true<-dplyr::filter(sum_cases_true, fx_etaria==a)
          nested_ordered_true<-nested_unordered_true[order(nested_unordered_true$onset_week),]
          # Sum nested counts over delay
          nested_unordered_true_sum<- nested_ordered_true|>dplyr::group_by(onset_week)|>dplyr::summarise(nested_col=sum(nested_col))
          # Fill in missing values:
          nested_list_true[[a]]<-tibble::tibble(onset_week=sort(rep(1:N_raw_true)))
          nested_list_true[[a]]<-dplyr::full_join(nested_list_true[[a]],nested_unordered_true_sum[,c("onset_week", "nested_col")],by = c("onset_week"))
          nested_list_true[[a]]$nested_col[is.na(nested_list_true[[a]]$nested_col)]<-0
          nested_list_true[[a]]<-nested_list_true[[a]]$nested_col
        }else{
          nested_unordered_true<-dplyr::filter(sum_cases_true, fx_etaria==a)
          nested_ordered_true<-nested_unordered_true[order(nested_unordered_true$onset_week),]
          # Sum nested counts over delay:
          nested_unordered_true_sum <- nested_ordered_true|>dplyr::group_by(onset_week,region_col)|>dplyr::summarise(nested_col=sum(nested_col))
          # Fill in missing values:
          nested_list_true[[a]]<-tibble::tibble(onset_week=sort(rep(1:N_raw_true,region_length_true)), region_col=rep(region_names,N_raw_true))
          nested_list_true[[a]]<-dplyr::full_join(nested_list_true[[a]],nested_unordered_true_sum[,c("onset_week", "nested_col","region_col")],by = c("onset_week", "region_col"))
          nested_list_true[[a]]$nested_col[is.na(nested_list_true[[a]]$nested_col)]<-0
          # Just nested count:
          nested_list_true[[a]]<-nested_list_true[[a]]$nested_col

        }

      }
    }


    # array final with dim by age:
    if(is.null(region_col)){
      if(is.null(age_col)){
        #array[TIME, DELAY]
        data.delay.array.available<-array(totals_list_true[[1]],dim=c(N_raw_true,(D_dim_true+1)))

        # Calculate observed totals
        data.totals.array.available<-apply(data.delay.array.available,c(1),sum)
        if(!is.null(nested_col)){
          data.nested.array.available<-nested_list_true[[1]]$nested_col
          print("Output data.nested.array.available has dimensions with the following labels data.nested.array.available[TIME].")

        }else{
          data.nested.array.available<-NULL
        }
        # Inform user of dimension order:
        print("Output data.delay.array.available has dimensions with the following labels data.delay.array.available[TIME, DELAY].")
        print("Output data.totals.array.available has dimensions with the following labels data.totals.array.available[TIME].")
        # Censor data.delay.array.available for run off triangle


      }else{
        #array[TIME, DELAY, AGE GROUP]
        data.delay.array.available<-array(abind::abind(totals_list_true[1:Age_length_true]),dim=c(N_raw_true,(D_dim_true+1),Age_length_true))

        # Calculate observed totals
        data.totals.array.available<-apply(data.delay.array.available,c(1,3),sum)
        # Inform user of dimension order:
        print("Output data.delay.array.available has dimensions with the following labels data.delay.array.available[TIME, DELAY, AGE GROUP].")
        print("Output data.totals.array.available has dimensions with the following labels data.totals.array.available[TIME, AGE GROUP].")
        if(!is.null(nested_col)){
          data.nested.array.available<-array(abind::abind(nested_list_true[1:Age_length_true]),dim=c(N_raw_true,Age_length_true))
          print("Output data.nested.array.available has dimensions with the following labels data.nested.array.available[TIME, AGE GROUP].")

        }else{
          data.nested.array.available<-NULL
        }

      }
    }else{
      if(is.null(age_col)){
        #array[REGION, TIME, DELAY]
        data.delay.array.available<-array(totals_list_true[[1]],dim=c(region_length_true, N_raw_true,(D_dim_true+1)))

        # Calculate observed totals
        #array[REGION, TIME]
        data.totals.array.available<-apply(data.delay.array.available,c(1,2),sum)

        # rearrange dimensions:
        #array[TIME, DELAY, REGION]
        data.delay.array.available<-aperm(data.delay.array.available, c(2,3,1))
        #array[TIME, REGION]
        data.totals.array.available<-aperm(data.totals.array.available, c(2,1))
        # Inform user of dimension order:
        print("Output data.delay.array.available has dimensions with the following labels data.delay.array.available[TIME, DELAY, REGION].")
        print("Output data.totals.array.available has dimensions with the following labels data.totals.array.available[TIME, REGION].")
        if(!is.null(nested_col)){
          data.nested.array.available<-array(abind::abind(nested_list_true[[1]]),dim=c(region_length_true, N_raw_true))
          data.nested.array.available<-aperm(data.nested.array.available, c(2,1))
          print("Output data.nested.array.available has dimensions with the following labels data.nested.array.available[TIME, REGION].")
        }else{
          data.nested.array.available<-NULL
        }



      }else{
        #array[REGION, TIME, DELAY, AGE GROUP]
        data.delay.array.available<-array(abind::abind(totals_list_true[1:Age_length_true]),dim=c(region_length_true, N_raw_true,(D_dim_true+1),Age_length_true))

        # Calculate observed totals
        data.totals.array.available<-apply(data.delay.array.available,c(1,2,4),sum)

        #array[TIME, DELAY, REGION, AGE GROUP]
        data.delay.array.available<-aperm(data.delay.array.available, c(2,3,1,4))
        #array[TIME, REGION, AGE GROUP]
        data.totals.array.available<-aperm(data.totals.array.available, c(2,1,3))
        # Inform user of dimension order:
        print("Output data.delay.array.available has dimensions with the following labels data.delay.array.available[TIME, DELAY, REGION, AGE GROUP].")
        print("Output data.totals.array.available has dimensions with the following labels data.totals.array.available[TIME, REGION, AGE GROUP].")
        if(!is.null(nested_col)){
          data.nested.array.available<-array(abind::abind(nested_list_true[1:Age_length_true]),dim=c(region_length_true,N_raw_true,Age_length_true))
          data.nested.array.available<-aperm(data.nested.array.available, c(2,1,3))
          print("Output data.nested.array.available has dimensions with the following labels data.nested.array.available[TIME, REGION, AGE GROUP].")
        }else{
          data.nested.array.available<-NULL
        }


      }
    }

  }else{
    data.totals.array.available<-NULL
    data.delay.array.available<-NULL
    data.nested.array.available<-NULL
  }


  ## Accounting for the maximum of days on the last week to be used
  data_w <- dataset |>
    dplyr::rename(date_report = {{date_report}},
                  date_onset = {{date_onset}},
                  age_col = {{age_col}},
                  nested_col ={{nested_col}},
                  region_col ={{region_col}}) |>
    dplyr::filter(date_report <= DT_max - aux.trimming.date,
                  date_onset>=(DT_max - aux.trimming.date - window.data.w + 1)) |>
    tidyr::drop_na(age_col) |>
    dplyr::mutate(
      ## Onset date
      # Moving the date to sunday
      DT.sun.aux = as.integer(format(date_onset, "%w")),
      ## Altering the date for the first day of the week
      dt.aux = date_onset -
        # Last recording date (DT_max_diadasemana) is the last day of the new week format
        DT.sun.aux +
        ifelse( use.epiweek, 0, DT_max_diadasemana+1 -
                  ifelse(DT_max_diadasemana+1>DT.sun.aux,7, 0)
        ),
      date_onset = dt.aux - ifelse( date_onset < dt.aux, 7, 0),
      # Recording date
      DT.sun.aux.rep = as.integer(format(date_report, "%w")),
      ## Altering the date for the first day of the week
      dt.aux.rep = date_report -
        # Last recording date (DT_max_diadasemana) is the last day of the new week format
        DT.sun.aux.rep +
        ifelse( use.epiweek, 0, DT_max_diadasemana+1 -
                  ifelse(DT_max_diadasemana+1 > DT.sun.aux.rep,7, 0)
        ),
      date_report = dt.aux.rep - ifelse( date_report < dt.aux.rep, 7, 0),
      Delay = as.numeric(date_report - date_onset) / 7,

    ) |>
    dplyr::select(-dt.aux,-dt.aux.rep, -DT.sun.aux,-DT.sun.aux.rep) |>
    dplyr::filter(Delay >= 0)


  ## Minimum date of onset
  first_date <- min(data_w |>
                      dplyr::pull(var = date_onset))

  ## Calculate onset weeks
  data_w<-dplyr::mutate(data_w, onset_week = (as.numeric((date_onset-first_date)/7)+1))

  # Create age col and remove any data greater than maximum age bin
  if(!is.null(age_col)){
    data_w <- data_w |>
      dplyr::mutate(age = floor(as.numeric(date_onset-age_col)/365.24),
             fx_etaria = cut(age,  breaks = bins_age, labels = labels_age, right = F))|>
      dplyr::filter(age <= max(bins_age)) |>
      tidyr::drop_na(fx_etaria)
    # Calculate number of cases for each week
    # Maximum regions
    if(!is.null(region_col)){
      sum_cases<-data_w|>dplyr::group_by(date_report,date_onset,onset_week,Delay,fx_etaria,region_col)|>dplyr::summarise(nested_col=sum(nested_col),cases=dplyr::n())
      sum_cases<-sum_cases|> dplyr::mutate(region_col=as.factor(region_col))|> dplyr::mutate(region=as.numeric(region_col))
      region_length<-max((sum_cases$region))
      region_names<-levels(sum_cases$region_col)
    }else{
      sum_cases<-data_w|>dplyr::group_by(date_report,date_onset,onset_week,Delay,fx_etaria)|>dplyr::summarise(nested_col=sum(nested_col),cases=dplyr::n())
      region_length=1
    }
    Age_length<-length(levels(sum_cases$fx_etaria))

  }else{
    # Calculate number of cases for each week
    if(!is.null(region_col)){
      sum_cases<-data_w|>dplyr::group_by(date_report,date_onset,onset_week,Delay,region_col)|>dplyr::summarise(nested_col=sum(nested_col),cases=dplyr::n())
      sum_cases<-sum_cases|> dplyr::mutate(region_col=as.factor(region_col))|> dplyr::mutate(region=as.numeric(region_col))
      region_length<-max((sum_cases$region))
      region_names<-levels(sum_cases$region_col)
    }else{
      sum_cases<-data_w|>dplyr::group_by(date_report,date_onset,onset_week,Delay)|>dplyr::summarise(nested_col=sum(nested_col),cases=dplyr::n())
      region_length=1
      region_names<-NULL
    }
    Age_length<-1
    sum_cases<-sum_cases|>dplyr::mutate(fx_etaria=Age_length)
  }
  if(!is.null(region_col)){
    sum_cases<-sum_cases[order(sum_cases$region_col),]

  }
  sum_cases<-sum_cases[order(sum_cases$onset_week),]

  # Maximum onset week
  N_raw <- max(sum_cases$onset_week)
  # Nowcast up to maximum reporting week
  N_raw <- N_raw + max(dplyr::filter(sum_cases,onset_week==max(sum_cases$onset_week))$Delay)

  # Maximum delay in data
  D_raw<-max(as.numeric(sum_cases$Delay),0)
  totals_list<-list()
  nested_list<-list()

  for(a in 1:Age_length){
    # Number of delays in data
    totals_list_unordered<-tidyr::pivot_wider(dplyr::select(dplyr::filter(sum_cases, fx_etaria==a),c(fx_etaria,date_report,date_onset,Delay,cases,onset_week,region_col)) ,id_cols=c(onset_week,region_col),names_from = Delay,values_from = cases)
    totals_list_unordered[is.na(totals_list_unordered)]<-0
    if(!is.null(region_col)){
      totals_list_unordered<-totals_list_unordered[order(totals_list_unordered$region_col),]

    }
    totals_list_unordered<-totals_list_unordered[order(totals_list_unordered$onset_week),]

    # Order delays
    order_delay<-totals_list_unordered[,-c(1:length(c(date_onset,region_col)))]
    order_delay<-order_delay[,order(as.integer(colnames(order_delay)))]
    totals_list_ordered1<-cbind(totals_list_unordered[,c(1:length(c(date_onset,region_col)))],order_delay)

    # Order attributes in rows:
    totals_list_ordered2<-totals_list_ordered1[,1:dim(totals_list_ordered1)[2]]|> dplyr::mutate(id=paste(onset_week,"-", region_col, sep=''))

    # Create blank data frame for all regions and weeks.
    if(is.null(region_col)){
      totals_list[[a]] <- cbind(tibble::tibble(onset_week=sort(rep((1:N_raw)))),
                                matrix(0, ncol=D_raw+1,nrow=N_raw))|>dplyr::mutate(check = NA)|>dplyr::mutate(id=paste(onset_week,sep=''))

    }else{
      totals_list[[a]] <- cbind(tibble::tibble(onset_week=sort(rep((1:N_raw),region_length)), region_col=rep(region_names,N_raw)),
                                matrix(0, ncol=D_raw+1,nrow=N_raw*region_length))|>dplyr::mutate(check = NA)|>dplyr::mutate(id=paste(onset_week,"-",region_col,sep=''))

    }

    # Fill in available data:
    totals_list[[a]][totals_list[[a]]$id %in% totals_list_ordered2$id, c(rep(FALSE,length(c(date_onset, region_col))),0:D_raw %in% as.numeric(colnames(order_delay)),TRUE,FALSE)] <- totals_list_ordered2[,(length(c(date_onset,region_col))+1):(dim(totals_list_ordered2)[2])]

    # Create D_max delays
    if(!is.null(D_max)){
      totals_list[[a]]<-cbind(totals_list[[a]][,(length(c(date_onset, region_col))+1):(D_max+length(c(date_onset, region_col)))],apply(totals_list[[a]][,(D_max+length(c(date_onset, region_col))+1):(dim(totals_list[[a]])[2]-2)],1,sum, na.rm=T))
      colnames(totals_list[[a]])[dim(totals_list[[a]])[2]]<-paste(D_max+1)
      D_dim<-D_max
    }else{
      D_dim<-D_raw
    }

    # Make a dimension for region if including region: array[REGION, TIME, DELAY]
    if(!is.null(region_col)){
      totals_list[[a]]<-totals_list[[a]][,1:(D_dim+1)]|>as.matrix()|>array(dim=c(region_length,N_raw,D_dim+1))
    }

  if(!is.null(nested_col)){
    if(is.null(region_col)){
      nested_unordered<-dplyr::filter(sum_cases, fx_etaria==a)
      nested_ordered<-nested_unordered[order(nested_unordered$onset_week),]
      # Sum nested counts over delay
      nested_unordered_sum<- nested_ordered|>dplyr::group_by(onset_week)|>dplyr::summarise(nested_col=sum(nested_col))
      # Fill in missing values:
      nested_list[[a]]<-tibble::tibble(onset_week=sort(rep(1:N_raw)))
      nested_list[[a]]<-dplyr::full_join(nested_list[[a]],nested_unordered_sum[,c("onset_week", "nested_col")],by = c("onset_week"))
      nested_list[[a]]$nested_col[is.na(nested_list[[a]]$nested_col)]<-0
       nested_list[[a]]<-nested_list[[a]]$nested_col
    }else{
      nested_unordered<-dplyr::filter(sum_cases, fx_etaria==a)
      nested_ordered<-nested_unordered[order(nested_unordered$onset_week),]
      # Sum nested counts over delay:
      nested_unordered_sum <- nested_ordered|>dplyr::group_by(onset_week,region_col)|>dplyr::summarise(nested_col=sum(nested_col))
      # Fill in missing values:
      nested_list[[a]]<-tibble::tibble(onset_week=sort(rep(1:N_raw,region_length)), region_col=rep(region_names,N_raw))
      nested_list[[a]]<-dplyr::full_join(nested_list[[a]],nested_unordered_sum[,c("onset_week", "nested_col","region_col")],by = c("onset_week", "region_col"))
      nested_list[[a]]$nested_col[is.na(nested_list[[a]]$nested_col)]<-0
      # Just nested count:
      nested_list[[a]]<-nested_list[[a]]$nested_col

    }

  }
  }


  # array final with dim by age:
  if(is.null(region_col)){
    if(is.null(age_col)){
      #array[TIME, DELAY]
      data.delay.array<-array(totals_list[[1]],dim=c(N_raw,(D_dim+1)))
      # Censor data.delay.array for run off triangle
      data.delay.array[outer(1:dim(data.delay.array[1:N_raw,])[1], 0:(dim(data.delay.array[1:N_raw,])[2]-1), FUN = "+") > N_raw] <- NA

      # Calculate observed totals
      data.totals.array<-apply(data.delay.array,c(1),sum)
      if(!is.null(nested_col)){
        data.nested.array<-nested_list[[1]]$nested_col
        print("Output data.nested.array has dimensions with the following labels data.nested.array[TIME].")

      }else{
        data.nested.array<-NULL
      }
      # Inform user of dimension order:
      print("Output data.delay.array has dimensions with the following labels data.delay.array[TIME, DELAY].")
      print("Output data.totals.array has dimensions with the following labels data.totals.array[TIME].")
      # Censor data.delay.array for run off triangle


      }else{
      #array[TIME, DELAY, AGE GROUP]
      data.delay.array<-array(abind::abind(totals_list[1:Age_length]),dim=c(N_raw,(D_dim+1),Age_length))
      # Censor data.delay.array for run off triangle
      for(a in 1:Age_length){
        data.delay.array[,,a][outer(1:dim(data.delay.array[1:N_raw,,a])[1], 0:(dim(data.delay.array[1:N_raw,,a])[2]-1), FUN = "+") > N_raw] <- NA
      }

    # Calculate observed totals
      data.totals.array<-apply(data.delay.array,c(1,3),sum)
      # Inform user of dimension order:
      print("Output data.delay.array has dimensions with the following labels data.delay.array[TIME, DELAY, AGE GROUP].")
      print("Output data.totals.array has dimensions with the following labels data.totals.array[TIME, AGE GROUP].")
      if(!is.null(nested_col)){
        data.nested.array<-array(abind::abind(nested_list[1:Age_length]),dim=c(N_raw,Age_length))
        print("Output data.nested.array has dimensions with the following labels data.nested.array[TIME, AGE GROUP].")

      }else{
        data.nested.array<-NULL
      }

    }
  }else{
    if(is.null(age_col)){
      #array[REGION, TIME, DELAY]
      data.delay.array<-array(totals_list[[1]],dim=c(region_length, N_raw,(D_dim+1)))
      # Censor data.delay.array for run off triangle
      for(s in 1:region_length){
        data.delay.array[s,,][outer(1:dim(data.delay.array[s,,])[1], 0:(dim(data.delay.array[s,,])[2]-1), FUN = "+") > N_raw] <- NA
      }
      # Calculate observed totals
      #array[REGION, TIME]
      data.totals.array<-apply(data.delay.array,c(1,2),sum)

      # rearrange dimensions:
      #array[TIME, DELAY, REGION]
      data.delay.array<-aperm(data.delay.array, c(2,3,1))
      #array[TIME, REGION]
      data.totals.array<-aperm(data.totals.array, c(2,1))
      # Inform user of dimension order:
      print("Output data.delay.array has dimensions with the following labels data.delay.array[TIME, DELAY, REGION].")
      print("Output data.totals.array has dimensions with the following labels data.totals.array[TIME, REGION].")
      if(!is.null(nested_col)){
        data.nested.array<-array(abind::abind(nested_list[[1]]),dim=c(region_length, N_raw))
        data.nested.array<-aperm(data.nested.array, c(2,1))
        print("Output data.nested.array has dimensions with the following labels data.nested.array[TIME, REGION].")
      }else{
        data.nested.array<-NULL
      }



    }else{
      #array[REGION, TIME, DELAY, AGE GROUP]
      data.delay.array<-array(abind::abind(totals_list[1:Age_length]),dim=c(region_length, N_raw,(D_dim+1),Age_length))
      # Censor data.delay.array for run off triangle
      for(s in 1:region_length){
        for(a in 1:Age_length){
          data.delay.array[s,,,a][outer(1:dim(data.delay.array[s,,,a])[1], 0:(dim(data.delay.array[s,,,a])[2]-1), FUN = "+") > N_raw] <- NA
        }
      }

      # Calculate observed totals
      data.totals.array<-apply(data.delay.array,c(1,2,4),sum)

       #array[TIME, DELAY, REGION, AGE GROUP]
      data.delay.array<-aperm(data.delay.array, c(2,3,1,4))
      #array[TIME, REGION, AGE GROUP]
      data.totals.array<-aperm(data.totals.array, c(2,1,3))
      # Inform user of dimension order:
      print("Output data.delay.array has dimensions with the following labels data.delay.array[TIME, DELAY, REGION, AGE GROUP].")
      print("Output data.totals.array has dimensions with the following labels data.totals.array[TIME, REGION, AGE GROUP].")
      if(!is.null(nested_col)){
        data.nested.array<-array(abind::abind(nested_list[1:Age_length]),dim=c(region_length,N_raw,Age_length))
        data.nested.array<-aperm(data.nested.array, c(2,1,3))
        print("Output data.nested.array has dimensions with the following labels data.nested.array[TIME, REGION, AGE GROUP].")
      }else{
        data.nested.array<-NULL
      }


    }
  }


  if(forecast>0){
    # Add NA's to delay matrix
    forecast.delay<-array(NA, dim=c(forecast,dim(data.delay.array)[-1]))
    data.delay.array<-abind::abind(data.delay.array,forecast.delay, along=1)
    # Add NA's to totals matrix
    forecast.totals<-array(NA, dim=c(forecast,dim(data.totals.array)[-1]))
    data.totals.array<-abind::abind(data.totals.array,forecast.totals, along=1)
    if(!is.null(nested_col)){
      data.nested.array<-abind::abind(data.nested.array,forecast.totals, along=1)
    }
    print(paste0("Setting forecast date as ", (first_date + (dim(data.totals.array)[1])*7),"."))
  }
  date_names<-seq.Date(first_date, (first_date + (dim(data.totals.array)[1])*7), by = "week")

  # Output
  return(list(data.totals.array=data.totals.array, data.delay.array=data.delay.array, data.nested.array=data.nested.array,
              data.totals.array.available=data.totals.array.available, data.delay.array.available=data.delay.array.available, data.nested.array.available=data.nested.array.available,
              region_names=region_names, date_names=date_names, bins_age=bins_age ))

}

