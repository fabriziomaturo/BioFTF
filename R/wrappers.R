## ---------------------------------------------------------------------------
## 7. Wrappers
## ---------------------------------------------------------------------------

export_plots <- function(x, file, from = 0, to = 2, n = 101) {
  if (missing(file) || !is.character(file) || length(file) != 1L) {
    stop("'file' must be the path of the PDF file to create.",
         call. = FALSE)
  }
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  hill_profile_plot(x, from = from, to = to, n = n)
  hill_derivative_plot(x, order = 1, from = from, to = to, n = n)
  hill_derivative_plot(x, order = 2, from = from, to = to, n = n)
  hill_curvature_plot(x, from = from, to = to, n = n)
  hill_radius_plot(x, from = from, to = to, n = n)
  hill_cum_plot(x, from = from, to = to, n = n)
  invisible(file)
}

alltools <- function(x, from = 0, to = 2, n = 101) {
  list(
    relative_abundances = datirel(x),
    species_summary = summary_species(x),
    hill_profile = hill_profile(x, from = from, to = to, n = n),
    first_derivative = hill_first_derivative(x, from = from, to = to, n = n),
    second_derivative = hill_second_derivative(x, from = from, to = to,
                                               n = n),
    curvature = hill_curvature(x, from = from, to = to, n = n),
    radius = hill_radius(x, from = from, to = to, n = n),
    area = hill_area(x, from = from, to = to, n = n),
    cumulative_integral = hill_cum(x, from = from, to = to, n = n),
    arc = hill_arc(x, from = from, to = to, n = n),
    ranking = hill_ranking(x, from = from, to = to, n = n)
  )
}

