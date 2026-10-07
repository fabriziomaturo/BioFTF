## ---------------------------------------------------------------------------
## 6. Biodiversity surfaces
## ---------------------------------------------------------------------------

.surface_time <- function(time, labels) {
  if (is.null(time)) {
    time <- suppressWarnings(as.numeric(labels))
    if (anyNA(time) || any(diff(time) <= 0)) {
      time <- seq_along(labels)
    }
  }
  if (!is.numeric(time) || length(time) != length(labels) ||
      any(!is.finite(time))) {
    stop("'time' must be a finite numeric vector with one value for each ",
         "row of the data.", call. = FALSE)
  }
  if (length(time) < 2L || any(diff(time) <= 0)) {
    stop("A surface requires at least two rows and strictly increasing ",
         "'time' values.", call. = FALSE)
  }
  as.numeric(time)
}

hill_surface <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                         time = NULL, allow_negative = FALSE) {
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  q <- attr(y, "q")
  .finite_grid(q, "Surfaces")
  time <- .surface_time(time, colnames(y))
  z <- matrix(as.vector(y), nrow = length(q),
              dimnames = list(.q_labels(q), colnames(y)))
  list(q = q, time = time, z = z)
}

.persp_surface <- function(s, xlab, ylab, zlab, main, theta, phi,
                           palette = "viridis", zlim = range(s$z),
                           diverging = FALSE, ...) {
  nq <- length(s$q)
  nt <- length(s$time)
  facet <- (s$z[-1L, -1L] + s$z[-1L, -nt] + s$z[-nq, -1L] +
              s$z[-nq, -nt]) / 4
  if (diff(zlim) <= 0) {
    zlim <- zlim + c(-0.5, 0.5)
  }
  crange <- if (diverging) c(-1, 1) * max(abs(zlim)) else zlim
  cols <- grDevices::hcl.colors(100L, palette)
  index <- 1L + floor(99 * (facet - crange[1L]) / diff(crange))
  index <- pmin(pmax(index, 1L), 100L)
  pmat <- graphics::persp(
    x = s$q, y = s$time, z = s$z, zlim = zlim, col = cols[index],
    border = NA, theta = theta, phi = phi, ticktype = "detailed",
    cex.axis = 0.8, xlab = paste0("\n", xlab), ylab = paste0("\n", ylab),
    zlab = paste0("\n", zlab), main = main, ...
  )
  for (j in seq_len(nt)) {
    graphics::lines(grDevices::trans3d(s$q, s$time[j], s$z[, j], pmat),
                    col = "grey20", lwd = 0.6)
  }
  for (i in unique(round(seq(1, nq, length.out = min(nq, 11L))))) {
    graphics::lines(grDevices::trans3d(s$q[i], s$time, s$z[i, ], pmat),
                    col = "grey20", lwd = 0.6)
  }
  invisible(pmat)
}

.need_plotly <- function() {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Interactive surfaces require the 'plotly' package. Install it ",
         "or use plot = \"static\".", call. = FALSE)
  }
}

.scene <- function(xlab, ylab, zlab) {
  list(xaxis = list(title = xlab), yaxis = list(title = ylab),
       zaxis = list(title = zlab))
}

hill_surface_plot <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                              time = NULL,
                              plot = c("interactive", "static"),
                              allow_negative = FALSE,
                              xlab = "Order q", ylab = "Time",
                              zlab = "Hill number",
                              main = "Hill-number surface",
                              theta = 40, phi = 25, ...) {
  plot <- match.arg(plot)
  s <- hill_surface(x, from = from, to = to, n = n, domain = domain,
                    time = time, allow_negative = allow_negative)
  if (plot == "interactive") {
    .need_plotly()
    p <- plotly::plot_ly(x = s$q, y = s$time, z = t(s$z), type = "surface",
                         colorscale = "Viridis", ...)
    return(plotly::layout(p, title = main,
                          scene = .scene(xlab, ylab, zlab)))
  }
  .persp_surface(s, xlab = xlab, ylab = ylab, zlab = zlab, main = main,
                 theta = theta, phi = phi, ...)
  invisible(s)
}

hill_surface_compare <- function(x, y, from = 0, to = 2, n = 101,
                                 domain = NULL, time = NULL,
                                 labels = c("Surface 1", "Surface 2"),
                                 plot = c("interactive", "static"),
                                 mode = c("overlay", "difference"),
                                 allow_negative = FALSE,
                                 main = "Comparison of Hill-number surfaces",
                                 theta = 40, phi = 25, ...) {
  plot <- match.arg(plot)
  mode <- match.arg(mode)
  if (length(labels) != 2L) {
    stop("'labels' must have length 2.", call. = FALSE)
  }
  if (NROW(x) != NROW(y)) {
    stop("'x' and 'y' must have the same number of rows.", call. = FALSE)
  }
  sx <- hill_surface(x, from = from, to = to, n = n, domain = domain,
                     time = time, allow_negative = allow_negative)
  sy <- hill_surface(y, from = from, to = to, n = n, domain = domain,
                     time = sx$time, allow_negative = allow_negative)
  sd <- list(q = sx$q, time = sx$time, z = sy$z - sx$z)
  out <- list(surface1 = sx, surface2 = sy, difference = sd)
  dmain <- paste0(main, ": ", labels[2L], " minus ", labels[1L])
  if (plot == "interactive") {
    .need_plotly()
    if (mode == "difference") {
      lim <- max(abs(sd$z))
      if (lim <= 0) lim <- 1
      p <- plotly::plot_ly(x = sd$q, y = sd$time, z = t(sd$z),
                           type = "surface", colorscale = "RdBu",
                           cmin = -lim, cmax = lim, ...)
      return(plotly::layout(p, title = dmain,
                            scene = .scene("Order q", "Time", "Difference")))
    }
    flat <- function(colour) list(list(0, colour), list(1, colour))
    p <- plotly::plot_ly(...)
    p <- plotly::add_surface(p, x = sx$q, y = sx$time, z = t(sx$z),
                             name = labels[1L], colorscale = flat("#1F77B4"),
                             opacity = 0.75, showscale = FALSE,
                             showlegend = TRUE)
    p <- plotly::add_surface(p, x = sy$q, y = sy$time, z = t(sy$z),
                             name = labels[2L], colorscale = flat("#D95F02"),
                             opacity = 0.75, showscale = FALSE,
                             showlegend = TRUE)
    return(plotly::layout(p, title = main,
                          scene = .scene("Order q", "Time", "Hill number")))
  }
  if (mode == "difference") {
    .persp_surface(sd, xlab = "Order q", ylab = "Time", zlab = "Difference",
                   main = dmain, theta = theta, phi = phi,
                   palette = "Blue-Red 3", diverging = TRUE, ...)
  } else {
    oldpar <- graphics::par(mfrow = c(1, 2), mar = c(2, 1, 3, 1))
    on.exit(graphics::par(oldpar), add = TRUE)
    zlim <- range(sx$z, sy$z)
    .persp_surface(sx, xlab = "Order q", ylab = "Time", zlab = "Hill number",
                   main = labels[1L], theta = theta, phi = phi,
                   zlim = zlim, ...)
    .persp_surface(sy, xlab = "Order q", ylab = "Time", zlab = "Hill number",
                   main = labels[2L], theta = theta, phi = phi,
                   zlim = zlim, ...)
  }
  invisible(out)
}
