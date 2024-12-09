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
#' @param date_format The format of the dates in data.individual.
#' @param forecast The number of weeks of forecasting you want to include (and to be passed to nimbleCast).
#' @param use.epiweek
#' @param K
#' @param trim.data
#'
#' @return
#' @export
#'
#' @examples
data_to_matrix <- function(data.individual,
                           window = NULL,
                           D_max= NULL,  # output matrix will have D_max+1 delays where D_max+1 is all cases reported after delay D_max
                           nested_col = c(FALSE, "default.covid", nested),
                           bins_age = bins_age,
                           date_onset,
                           date_report,
                           age_col=NULL, # Allow for age as well as DOB?
                           region_col=NULL,
                           date_format=c("%d/%m/%Y"),
                           use.epiweek = FALSE,
                           K = 0,
                           trim.data){

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
    dataset<-dataset|>
      dplyr::select(-PCR_SARS2,-AN_SARS2)
  }

  # Set bins for ages:
  if(!is.null(age_col)){
    if(!is.numeric(bins_age)){
      if(missing(bins_age) | is.null(bins_age)){
        bins_age<-c(0,5,12,18,30,seq(40,90,by=10),130)
        labels_age<-1:(length(bins_age)-1)
        warning("Bins age as in SI-PNI: ",
                stringr::str_c(bins_age[bins_age < 90], " "), "90+",
                call. = T)
      } else if(bins_age == "SI-PNI"){
        bins_age<-c(0,5,12,18,30,seq(40,90,by=10),130)
        labels_age<-1:(length(bins_age)-1)
        warning("Bins age as in SI-PNI: ",
                stringr::str_c(bins_age[bins_age < 90], " "), "90+",
                call. = T)
      } else if(bins_age == "10 years"){
        bins_age<-c(seq(0,90, by = 10),130)
        labels_age<-1:(length(bins_age)-1)
        warning("Bins age in 10 years: ",
                stringr::str_c(bins_age[bins_age < 90], " "), "90+",
                call. = T)
      } else if(bins_age == "5 years"){
        bins_age<-c(seq(0,90, by = 5),130)
        labels_age<-1:(length(bins_age)-1)
        warning("Bins age in 5 years: ",
                stringr::str_c(bins_age[bins_age < 90], " "), "90+",
                call. = T)
      }
      else {
        stop("Age bins are not options of 'SI-PNI', '10 years', '5 years' or numeric vector!")
      }
    } else {
      bins_age<-bins_age
      labels_age<-1:(length(bins_age)-1)
      warning("Using bins ages given: ",
              stringr::str_c(bins_age[bins_age < bins_age[length(bins_age) - 1]], " "),
              bins_age[length(bins_age) - 1], "+",
              call. = T)
    }
  }else{
    bins_age<-c(0,130)
    labels_age<-1
    warning("No age_col given so data will be returned with no age groups")
  }


  ## Last digitation date considered
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



  ## K parameter of King
  K.w<-7*K

  ## Maximum date to be considered on the estimation
  DT_max <- max(dataset |>
                  dplyr::pull(var = {{date_report}}),
                na.rm = T) - trim.data.w + K.w
  ## Maximum date of onset
  DT_max_onset <- max(dataset |>
                        dplyr::pull(var = {{date_onset}}))

  ## Transforming window into days
  if(is.null(window)){
    window.data.w<- max(dataset |>
                          dplyr::pull(var = {{date_onset}}),
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


  ## Accounting for the maximum of days on the last week to be used
  data_w <- dataset |>
    dplyr::rename(date_report = {{date_report}},
                  date_onset = {{date_onset}},
                  age_col = {{age_col}},
                  nested_col ={{nested_col}},
                  region_col ={{region_col}}) |>
    dplyr::filter(date_report <= DT_max - aux.trimming.date,
                  date_onset>=DT_max_onset-window.data.w) |>
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
  data_w<-mutate(data_w, onset_week = (as.numeric((date_onset-first_date)/7)+1))

  # Create age col and remove any data greater than maximum age bin
  if(!is.null(age_col)){
    data_w <- data_w |>
      mutate(age = floor(as.numeric(date_onset-age_col)/365.24),
             fx_etaria = cut(age,  breaks = bins_age, labels = labels_age, right = F))|>
      filter(age <= max(bins_age)) |>
      tidyr::drop_na(fx_etaria)
    # Calculate number of cases for each week
    # Maximum regions
    if(!is.null(region_col)){
      sum_cases<-data_w%>%group_by(date_report,date_onset,onset_week,Delay,fx_etaria,region_col)%>%summarise(nested_col=sum(nested_col),cases=n())
      sum_cases<-sum_cases|> mutate(region_col=as.factor(region_col))|> mutate(region=as.numeric(region_col))
      region_length<-max((sum_cases$region))
      region_names<-levels(sum_cases$region_col)
    }else{
      sum_cases<-data_w%>%group_by(date_report,date_onset,onset_week,Delay,fx_etaria)%>%summarise(nested_col=sum(nested_col),cases=n())
      region_length=1
    }
    Age_length<-length(levels(sum_cases$fx_etaria))

  }else{
    # Calculate number of cases for each week
    if(!is.null(region_col)){
      sum_cases<-data_w%>%group_by(date_report,date_onset,onset_week,Delay,region_col)%>%summarise(nested_col=sum(nested_col),cases=n())
      sum_cases<-sum_cases|> mutate(region_col=as.factor(region_col))|> mutate(region=as.numeric(region_col))
      region_length<-max((sum_cases$region))
      region_names<-levels(sum_cases$region_col)
    }else{
      sum_cases<-data_w%>%group_by(date_report,date_onset,onset_week,Delay)%>%summarise(nested_col=sum(nested_col),cases=n())
      sum_cases<-data_w%>%group_by(date_report,date_onset,onset_week,Delay,fx_etaria)%>%summarise(nested_col=sum(nested_col),cases=n())
      region_length=1
    }
    Age_length<-1
    sum_cases<-sum_cases|>mutate(fx_etaria=Age_length)
  }

  # Maximum time steps in data
  N_raw<-max(sum_cases$onset_week)

  # Maximum delay in data
  D_raw<-max(as.numeric(sum_cases$Delay),0)
  totals_list<-list()
  nested_list<-list()
  for(a in 1:Age_length){
    # Number of delays in data
    totals_list_unordered<-pivot_wider(select(filter(sum_cases, fx_etaria==a),c(fx_etaria,date_report,date_onset,Delay,cases,onset_week,region_col)) ,id_cols=c(onset_week,region_col),names_from = Delay,values_from = cases)
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
    totals_list_ordered2<-totals_list_ordered1[,1:dim(totals_list_ordered1)[2]]%>%mutate(id=paste(onset_week,"-", region_col, sep=''))

    # Create blank data frame for all regions and weeks.
    if(is.null(region_col)){
      totals_list[[a]] <- cbind(tibble(onset_week=sort(rep((1:N_raw)))),
                                matrix(0, ncol=D_raw+1,nrow=N_raw))%>%mutate(check = NA)%>%mutate(id=paste(onset_week,sep=''))

    }else{
      totals_list[[a]] <- cbind(tibble(onset_week=sort(rep((1:N_raw),region_length)), region_col=rep(region_names,N_raw)),
                                matrix(0, ncol=D_raw+1,nrow=N_raw*region_length))%>%mutate(check = NA)%>%mutate(id=paste(onset_week,"-",region_col,sep=''))

    }

    # Fill in available data:
    totals_list[[a]][totals_list[[a]]$id %in% totals_list_ordered2$id, c(rep(FALSE,length(c(date_onset, region_col))),0:D_raw %in% as.numeric(colnames(order_delay)),TRUE,FALSE)] <- totals_list_ordered2[,(length(c(date_onset,region_col))+1):(dim(totals_list_ordered2)[2])]

    # Create D_max delays
    if(!is.null(D_max)){
      totals_list[[a]]<-cbind(totals_list[[a]][,(length(c(date_onset, region_col))+1):(D_max+length(c(date_onset, region_col)))],apply(totals_list[[a]][,(D_max+length(c(date_onset, region_col))+1):(dim(totals_list[[a]])[2]-2)],1,sum))
      colnames(totals_list[[a]])[dim(totals_list[[a]])[2]]<-paste(D_max+1)
      D_dim<-D_max
    }else{
      D_dim<-D_raw
    }

    # Make a dimension for region if including region: array[REGION, TIME, DELAY]
    if(!is.null(region_col)){
      totals_list[[a]]<-totals_list[[a]][,1:(D_dim+1)]%>%as.matrix()%>%array(dim=c(region_length,N_raw,D_dim+1))
    }
  }

  if(!is.null(nested_col)){
    if(is.null(region_col)){
      nested_unordered<-filter(sum_cases, fx_etaria==a)
      nested_ordered<-nested_unordered[order(nested_unordered$onset_week),]
      nested_list[[a]]<-tibble(onset_week=sort(rep(1:N_raw)))
      nested_list[[a]]<-full_join(nested_list[[a]],nested_ordered[,c("onset_week", "nested_col","Delay")],by = c("onset_week"))
      nested_list[[a]]$nested_col[is.na(nested_list[[a]]$nested_col)]<-0
      # Sum nested counts over delay
      nested_list[[a]]<- nested_list[[a]]%>%group_by(onset_week)%>%summarise(nested_col=sum(nested_col))
      nested_list[[a]]<-nested_list[[a]]$nested_col
    }else{
      nested_unordered<-filter(sum_cases, fx_etaria==a)
      nested_ordered<-nested_unordered[order(nested_unordered$onset_week),]
      nested_list[[a]]<-tibble(onset_week=sort(rep(1:N_raw,region_length)), region=rep(1:region_length,N_raw))
      nested_list[[a]]<-full_join(nested_list[[a]],nested_ordered[,c("onset_week", "region", "nested_col","region_col","Delay")],by = c("onset_week", "region"))
      nested_list[[a]]$nested_col[is.na(nested_list[[a]]$nested_col)]<-0
      # Sum nested counts over delay
      nested_list[[a]]<- nested_list[[a]]%>%group_by(onset_week,region)%>%summarise(nested_col=sum(nested_col))
      nested_list[[a]]<-nested_list[[a]]$nested_col
    }

  }


  # array final with dim by age:
  if(is.null(region_col)){
    if(is.null(age_col)){
      #array[TIME, DELAY]
      data.delay.array<-array(totals_list[[1]],dim=c(N_raw,(D_dim+1)))
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
    }else{
      #array[TIME, DELAY, AGE GROUP]
      data.delay.array<-array(abind(totals_list[1:Age_length]),dim=c(N_raw,(D_dim+1),Age_length))
      data.totals.array<-apply(data.delay.array,c(1,3),sum)
      # Inform user of dimension order:
      print("Output data.delay.array has dimensions with the following labels data.delay.array[TIME, DELAY, AGE GROUP].")
      print("Output data.totals.array has dimensions with the following labels data.totals.array[TIME, AGE GROUP].")
      if(!is.null(nested_col)){
        data.nested.array<-array(abind(nested_list[1:Age_length]),dim=c(N_raw,Age_length))
        print("Output data.nested.array has dimensions with the following labels data.nested.array[TIME, AGE GROUP].")

      }else{
        data.nested.array<-NULL
      }
    }

  }else{
    if(is.null(age_col)){
      #array[REGION, TIME, DELAY]
      data.delay.array<-array(totals_list[[1]],dim=c(region_length, N_raw,(D_dim+1)))
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
        data.nested.array<-array(abind(nested_list[[1]]),dim=c(region_length, N_raw))
        data.nested.array<-aperm(data.nested.array, c(2,1))
        print("Output data.nested.array has dimensions with the following labels data.nested.array[TIME, REGION].")

      }else{
        data.nested.array<-NULL
      }
    }else{
      #array[REGION, TIME, DELAY, AGE GROUP]
      data.delay.array<-array(abind(totals_list[1:Age_length]),dim=c(region_length, N_raw,(D_dim+1),Age_length))
      data.totals.array<-apply(data.delay.array,c(1,2,4),sum)
      #array[TIME, DELAY, REGION, AGE GROUP]
      data.delay.array<-aperm(data.delay.array, c(2,3,1,4))
      #array[TIME, REGION, AGE GROUP]
      data.totals.array<-aperm(data.totals.array, c(2,1,3))
      # Inform user of dimension order:
      print("Output data.delay.array has dimensions with the following labels data.delay.array[TIME, DELAY, REGION, AGE GROUP].")
      print("Output data.totals.array has dimensions with the following labels data.totals.array[TIME, REGION, AGE GROUP].")
      if(!is.null(nested_col)){
        data.nested.array<-array(abind(nested_list[1:Age_length]),dim=c(region_length,N_raw,Age_length))
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
    data.delay.array<-abind(data.delay.array,forecast.delay, along=1)
    # Add NA's to totals matrix
    forecast.totals<-array(NA, dim=c(forecast,dim(data.totals.array)[-1]))
    data.totals.array<-abind(data.totals.array,forecast.totals, along=1)
    if(!is.null(nested_col)){
      data.nested.array<-abind(data.nested.array,forecast.totals, along=1)
    }

  }


  # Output
  return(list(data.totals.array=data.totals.array, data.delay.array=data.delay.array, data.nested.array=data.nested.array))

}



