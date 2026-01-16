#' Plot the On/Off Kinetics of the provided data set
#'
#' @param rawdata rawdata object
#' @param type type of data to plot. Default setting will facet the plot
#' by type.
#'
#' @return ggplot
#' @export
#'
plot_rawdata <- function(rawdata, type = NULL) {
  # TODO stop if not object of class rawdata

  plotdata <- prepare_for_ook_plot(rawdata, type = NULL)
  # TODO add title?
  ggook(plotdata, x = NULL, y = NULL)
}

#' Manipulate data frame for use in plotting functions
#'
#' @param rawdata rawdata object
#' @param type type of data to display
#'
#' internal function for use in ook_plot and plot_rawdata
prepare_for_ook_plot <- function(rawdata, type = NULL) {
  prepare_shiny_plot(
    rawdata$data,
    rawdata$info$sfreq,
    rawdata$info$bounds$meas_start,
    rawdata$info$bounds$meas_end,
    rawdata$info$bads,
    rawdata$info$bounds$removed,
    rawdata$info$cols$col_type,
    rawdata$info$cols$col_name,
    type
  )
}

#' Wrapper around ggplot functions to plot on/off kinetics
#'
#' @param data data frame to use
#' @param x limits
#' @param y limits
ggook <- function(data, x, y) {
  cols <- c("ZeroedTime", "value", "col_name", "col_type")

  if (length(intersect(cols, colnames(data))) != length(cols)) {
    stop("Provided data.frame is not compatible with this function.")
  }
  major_breaks <- seq(0, 960, 120)
  minor_breaks <- seq(0, 960, 60)

  # TODO make manual color scale work
  plt <- ggplot(
    data,
    aes(x = ZeroedTime, y = value)
  ) +
    geom_line(aes(color = col_name), na.rm = TRUE) +
    theme_bw() +
    coord_cartesian(xlim = x, ylim = y, expand = FALSE) +
    labs(
      x = "Time (sec)",
      y = "",
      color = ""
    )

  if (is.null(x) & is.null(y)) {
    return(plt + scale_x_continuous(breaks = seq(0, 960, 120)))
  }

  bk <- major_breaks[which(between(major_breaks, x[1], x[2]))]
  vlines <- lapply(bk, function(x) geom_vline(xintercept = x))

  plt <- lapply(vlines, function(x) plt + x)

  plt
}

#' manual color scale for ggplot
#'
#' @param col_names character vector of unique column names in the dataset
#'
#' @return named vector of colors to use with each col_name
manual_color <- function(col_names) {
  manual_color_opts <- c(
    "S1_D1_759" = "firebrick2",
    "S1_D1_854" = "royalblue1",
    "S2_D1_759" = "red2",
    "S2_D1_854" = "blue",
    "S3_D1_761" = "red4",
    "S3_D1_856" = "royalblue4",
    "TSI" = "springgreen2",
    "TSIFF" = "springgreen4"
  )

  manual_color_opts[col_names]
}

#' Plot Regressions on top of actual values
#'
#' @param models analyzed object
#'
#' @return ggplot object
#' @export
plot_regressions <- function(models) {
  models$data |>
    dplyr::filter(Section != "BegRest" & Section != "WarmUp") |>
    ggplot(aes(x = SectionZeroedTime)) +
    facet_grid(col_type ~ Section, scales = "free_y") +
    geom_line(aes(y = value_norm), color = "grey23", linewidth = 1.25, alpha = 0.5) +
    geom_line(aes(y = predicted_norm), color = "red", linewidth = 1.25, na.rm = TRUE) +
    labs(
      x = "Time Since Start of Section (seconds)",
      y = ""
    ) +
    scale_x_continuous(breaks = c(0, 30, 60, 90, 120), limits = c(0, 120)) +
    scale_color_manual(
      labels = c("Actual", "Regression"),
      values = c("grey23", "red")
    ) +
    theme_light()
}
