#' nls wrapped with possibly
#'
#' @return function
#' @export
#' @noRd
possibly_nls <- purrr::possibly(nls, otherwise = NA)


#' Self Starting Asymptotic Regression
#'
#' @param df data frame
#'
#' @return list of coefficients from regression if successful, NA if unable to run
#' @export
asym_model <- function(df) {
  ## TODO use tidy evaluation and data masking to use input for x and y variables
  suppressWarnings(
    possibly_nls(value ~ SSasymp(SectionZeroedTime, Asym, R0, lrc), df)
  )
}

#' Calculate Normalized Asymptotic Values
#'
#' Subtracts R0 (starting value) from the predicted curve.
#' y_norm = y(t) - R0
#'
#' @param t time vector
#' @param Asym Asymptote
#' @param R0 Starting value
#' @param lrc log rate constant
#'
#' @export
predict_asymp_normalized <- function(t, Asym, R0, lrc) {
  # y(t) = Asym + (R0 - Asym) * exp(-exp(lrc) * t)
  # y(t) - R0
  y <- Asym + (R0 - Asym) * exp(-exp(lrc) * t)
  return(y - R0)
}


#' Run Regressions on On/Off Kinetics Data
#'
#' @param rawdata rawdata object
#'
#' @return dataframe
#' @export
run_models <- function(rawdata) {
  nested <- group_nest_dt(Section, col_type) |>
    setkey(Section, col_type)

  nested[, Model := map(data, asym_model)]
}

#' Calculate steady state values
#'
#' Default width used to calculate steady state is 30 seconds
#'
#' @param dt data.frame returned from `prepare_for_modeling()`
#' @param width A double, number of seconds
#'
#' @return data.table
#' @export
steady_states <- function(dt, width = 30) {
  dt[
    Section != "BegRest" & Section != "WarmUp" & SectionZeroedTime > val,
    .(steadystate = mean(value, na.rm = TRUE)),
    by = .(Section, col_type),
    env = list(val = width)
  ]
}
