% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_to_matrix.R
\name{data_to_matrix}
\alias{data_to_matrix}
\title{Convert individual patient data to a matrix.}
\usage{
data_to_matrix(
  data.individual,
  window = NULL,
  D_max = NULL,
  nested_col = c(FALSE, "default.covid", nested_col),
  bins_age = c(NULL, "SI-PNI", "10 years", "5 years", bins_age),
  date_onset,
  date_report,
  age_col = NULL,
  region_col = NULL,
  date_format = c("\%d/\%m/\%Y"),
  forecast = 0,
  use.epiweek = FALSE,
  trim.data,
  testing = FALSE
)
}
\arguments{
\item{data.individual}{Data set formatted with individual patient level data.}

\item{window}{Number of weeks from the most recent date you want to keep.}

\item{D_max}{Number of delays in weeks you want to be returned. Output will have D_max+1 delay dimensions where D_max+1 is all remainder counts after D_max.}

\item{nested_col}{Name of the column related to a potential nested structure you wish to model. Can be set to "default.covid" if COVID-19 test results are in columns "AN_SARS2" & "PCR_SARS2".}

\item{bins_age}{A numeric vector of bins you wish to divide age into e.g. c(0,25,50,75,100,125), or can be set to "SI-PNI", "5 years" or "10 year".}

\item{date_onset}{The column name corresponding to the date of onset.}

\item{date_report}{The column name corresponding to the date of reporting.}

\item{age_col}{The column name corresponding to the date of birth.}

\item{region_col}{The column name corresponding to the spatial region.}

\item{date_format}{The format of the dates in data.individual, default is "\%d/\%m/\%Y".}

\item{forecast}{Integer of weeks of forecasting you want to include (to be passed to nimbleCast).}

\item{use.epiweek}{If TRUE, it uses the CDC epiweek definition where the week starts on Sunday, if FALSE it the week ends at the weekday of the last record date.}

\item{trim.data}{Integer such that nowcasting is performed up to trim.data weeks prior to the maximum date in the data.}

\item{testing}{Default is FALSE which implys data is most recent data available and nowcast date is most recent date_report. Set to TRUE if you wish to nowcast up to the most recent date_onset. If set as date in format yyyy-mm-dd that will be the nowcast date. If not FALSE arrays with all available data will also be returned.}
}
\value{
List of data arrays.
}
\description{
Convert individual patient data to a matrix.
}
\examples{
Sari to matrix
}
