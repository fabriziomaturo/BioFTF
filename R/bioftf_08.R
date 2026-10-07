## ---------------------------------------------------------------------------
## 8. Deprecated names of BioFTF 1.x
## ---------------------------------------------------------------------------

.legacy <- function(old, new) {
  .Deprecated(new, package = "BioFTF", old = old, msg = paste0(
    "'", old, "()' is deprecated; use '", new, "()' instead. ",
    "Since BioFTF 3.0.0 this tool is computed on the Hill-number profile ",
    "(order q in [0, 2] by default) and no longer on the Patil-Taillie ",
    "beta profile (beta in [-1, 1]), so the values differ from those of ",
    "BioFTF 1.x."
  ))
}

beta <- function(x, n = 101, ...) {
  .legacy("beta", "hill_profile")
  hill_profile(x, n = n, ...)
}

beta_plot <- function(x, n = 101, ...) {
  .legacy("beta_plot", "hill_profile_plot")
  hill_profile_plot(x, n = n, ...)
}

betaprime <- function(x, n = 101, ...) {
  .legacy("betaprime", "hill_first_derivative")
  hill_first_derivative(x, n = n, ...)
}

betaprime_plot <- function(x, n = 101, ...) {
  .legacy("betaprime_plot", "hill_derivative_plot")
  hill_derivative_plot(x, order = 1, n = n, ...)
}

betasecond <- function(x, n = 101, ...) {
  .legacy("betasecond", "hill_second_derivative")
  hill_second_derivative(x, n = n, ...)
}

betasecond_plot <- function(x, n = 101, ...) {
  .legacy("betasecond_plot", "hill_derivative_plot")
  hill_derivative_plot(x, order = 2, n = n, ...)
}

curvature <- function(x, n = 101, ...) {
  .legacy("curvature", "hill_curvature")
  hill_curvature(x, n = n, ...)
}

curvature_plot <- function(x, n = 101, ...) {
  .legacy("curvature_plot", "hill_curvature_plot")
  hill_curvature_plot(x, n = n, ...)
}

radius <- function(x, n = 101, ...) {
  .legacy("radius", "hill_radius")
  hill_radius(x, n = n, ...)
}

radius_plot <- function(x, n = 101, ...) {
  .legacy("radius_plot", "hill_radius_plot")
  hill_radius_plot(x, n = n, ...)
}

arc <- function(x, n = 101, ...) {
  .legacy("arc", "hill_arc")
  hill_arc(x, n = n, ...)
}

area <- function(x, n = 101, ...) {
  .legacy("area", "hill_area")
  hill_area(x, n = n, ...)
}

ranking <- function(x, n = 101, ...) {
  .legacy("ranking", "hill_ranking")
  hill_ranking(x, n = n, ...)
}
