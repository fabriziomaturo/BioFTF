## ---------------------------------------------------------------------------
## 5. Profile plots
## ---------------------------------------------------------------------------

.curve_colours <- function(k) {
  if (k == 1L) "black" else grDevices::hcl.colors(k, "Dark 3")
}

.plot_curves <- function(y, xlab, ylab, main, legend = "topright", ...) {
  q <- attr(y, "q")
  if (any(!is.finite(q))) {
    stop("Plots require finite orders.", call. = FALSE)
  }
  k <- ncol(y)
  cols <- .curve_colours(k)
  ltys <- rep_len(1:6, k)
  graphics::matplot(q, y, type = "l", lty = ltys, col = cols, lwd = 2,
                    xlab = xlab, ylab = ylab, main = main, ...)
  if (!is.null(legend)) {
    graphics::legend(legend, legend = colnames(y), lty = ltys, col = cols,
                     lwd = 2, cex = 0.8, bty = "n")
  }
  invisible(y)
}

hill_profile_plot <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                              allow_negative = FALSE, legend = "topright",
                              ...) {
  .plot_curves(
    hill_profile(x, from = from, to = to, n = n, domain = domain,
                 allow_negative = allow_negative),
    xlab = "Order q", ylab = "Hill number",
    main = "Hill-number diversity profiles", legend = legend, ...
  )
}

hill_plot <- hill_profile_plot

hill_derivative_plot <- function(x, order = 1, from = 0, to = 2, n = 101,
                                 domain = NULL, allow_negative = FALSE,
                                 legend = "bottomright", ...) {
  y <- hill_derivative(x, order = order, from = from, to = to, n = n,
                       domain = domain, allow_negative = allow_negative)
  label <- if (order == 1) "First derivative" else "Second derivative"
  .plot_curves(y, xlab = "Order q", ylab = label,
               main = paste(label, "of Hill-number profiles"),
               legend = legend, ...)
}

hill_curvature_plot <- function(x, from = 0, to = 2, n = 101,
                                domain = NULL, allow_negative = FALSE,
                                legend = "topright", ...) {
  .plot_curves(
    hill_curvature(x, from = from, to = to, n = n, domain = domain,
                   allow_negative = allow_negative),
    xlab = "Order q", ylab = "Curvature",
    main = "Curvature of Hill-number profiles", legend = legend, ...
  )
}

hill_radius_plot <- function(x, from = 0, to = 2, n = 101,
                             domain = NULL, allow_negative = FALSE,
                             legend = "topleft", ...) {
  .plot_curves(
    hill_radius(x, from = from, to = to, n = n, domain = domain,
                allow_negative = allow_negative),
    xlab = "Order q", ylab = "Radius of curvature (log scale)",
    main = "Radius of curvature of Hill-number profiles",
    legend = legend, log = "y", ...
  )
}

hill_cum_plot <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                          baseline = 0, allow_negative = FALSE,
                          legend = "topleft", ...) {
  .plot_curves(
    hill_cum(x, from = from, to = to, n = n, domain = domain,
             baseline = baseline, allow_negative = allow_negative),
    xlab = "Order q", ylab = "Cumulative integral",
    main = "Cumulative integral functions of Hill-number profiles",
    legend = legend, ...
  )
}

