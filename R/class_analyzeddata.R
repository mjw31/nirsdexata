#' Create an Analyzed Data Object
#'
#' @param data Processed data frame with original values, predictions, and residuals
#' @param info Metadata list from the original rawdata object
#' @param coefs Data frame of regression coefficients (Asym, R0, lrc) and Tau
#' @param steadystates Data frame of steady state values by Section and col_type
#'
#' @return An object of class "analyzed"
#' @export
new_analyzed_data <- function(data,
                              info,
                              coefs = NULL,
                              steadystates = NULL) {
  stopifnot(is.data.frame(data))
  stopifnot(is.list(info))
  stopifnot(is.data.frame(coefs) || is.null(coefs))
  stopifnot(is.data.frame(steadystates) || is.null(steadystates))

  analyzed <- list(
    data = data,
    info = info,
    coefs = coefs,
    steadystates = steadystates
  )

  structure(analyzed, class = c("analyzed", "list"))
}

#' Analyzed Data helper
#'
#' @param rawdata An object of class "rawdata" containing data and metadata
#'
#' @return An object of class "analyzed" with model fits and steady states
#' @export
analysis <- function(rawdata) {
  prep <- prepare_for_modeling(rawdata) |> data.table::setDT()

  # Exclude initial sections before modeling
  prep2 <- prep[Section != "BegRest" & Section != "WarmUp", ]

  # Nest data by section and type
  nested <- group_nest_dt(prep2, Section, col_type) |> data.table::setkey(Section, col_type)

  # Fit asymptotic models and tidy results
  asym_models <- nested[, Model := purrr::map(data, asym_model)][
    !is.na(Model),
    `:=`(
      Tidy = purrr::map(Model, broom::tidy),
      Augment = purrr::map(Model, broom::augment)
    )
  ]

  # Remove any models that produced empty tidy outputs
  asym_models <- asym_models[lengths(Tidy) > 0]

  # Compute steady state values
  steadystates <- steady_states(prep2)

  # Unnest coefficients and merge back into nested structure
  tidy <- unnest_dt(asym_models,
    col = Tidy,
    id  = list(Section, col_type)
  ) |>
    data.table::dcast(Section + col_type ~ term,
      value.var = "estimate"
    ) |>
    merge(nested, all.y = TRUE)

  # Compute Tau and merge with steady states
  coefs <- tidy[, .(Section, col_type, Asym, R0, lrc)][, Tau := (1 / exp(lrc))] |>
    merge(y = steadystates)

  # Calculate steadystate_norm = steadystate - R0
  coefs[, steadystate_norm := steadystate - R0]

  # Update steadystates to include the normalized value
  steadystates <- coefs[, .(Section, col_type, steadystate, steadystate_norm)]

  # Prepare augmented predictions and residuals
  augment <- unnest_dt(asym_models,
    col = Augment,
    id  = list(Section, col_type)
  ) |>
    dplyr::rename(
      predicted = .fitted,
      residual = .resid
    ) |>
    dplyr::select(Section, col_type, SectionZeroedTime, predicted, residual)

  # Merge augmented data with original prepped data
  dt <- merge(prep,
    augment,
    by    = c("Section", "col_type", "SectionZeroedTime"),
    all.x = TRUE
  )

  # Join with coefficients to calculate normalized predictions
  dt <- merge(dt,
    coefs[, .(Section, col_type, Asym, R0, lrc)],
    by = c("Section", "col_type"),
    all.x = TRUE
  )

  # Calculate predicted_norm = predicted - R0
  # Note: predicted comes from nls, which is Asym + (R0-Asym)*exp(...).
  # So predicted - R0 is the normalized curve starting at 0.
  dt[, predicted_norm := predicted - R0]

  # Calculate value_norm = value - R0 for plotting normalized data
  dt[, value_norm := value - R0]

  # Return analyzed object including steady states
  new_analyzed_data(dt,
    rawdata$info,
    coefs        = coefs,
    steadystates = steadystates
  )
}
