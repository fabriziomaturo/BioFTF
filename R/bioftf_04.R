## ---------------------------------------------------------------------------
## 4. Profiles and functional tools
## ---------------------------------------------------------------------------

hill_profile <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                         allow_negative = FALSE) {
  q <- .hill_domain(domain = domain, from = from, to = to, n = n,
                    allow_negative = allow_negative)
  .hill_matrix(x, q, which = 1L)
}

hill <- hill_profile

hill_derivative <- function(x, order = 1, from = 0, to = 2, n = 101,
                            domain = NULL, allow_negative = FALSE) {
  if (length(order) != 1L || !order %in% c(1, 2)) {
    stop("'order' must be 1 or 2.", call. = FALSE)
  }
  q <- .hill_domain(domain = domain, from = from, to = to, n = n,
                    allow_negative = allow_negative)
  if (any(!is.finite(q))) {
    stop("Derivatives require finite orders.", call. = FALSE)
  }
  .hill_matrix(x, q, which = as.integer(order) + 1L)
}

hill_first_derivative <- function(x, from = 0, to = 2, n = 101,
                                  domain = NULL, allow_negative = FALSE) {
  hill_derivative(x, order = 1, from = from, to = to, n = n,
                  domain = domain, allow_negative = allow_negative)
}

hill_second_derivative <- function(x, from = 0, to = 2, n = 101,
                                   domain = NULL, allow_negative = FALSE) {
  hill_derivative(x, order = 2, from = from, to = to, n = n,
                  domain = domain, allow_negative = allow_negative)
}

hill_curvature <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                           allow_negative = FALSE) {
  d1 <- hill_derivative(x, order = 1, from = from, to = to, n = n,
                        domain = domain, allow_negative = allow_negative)
  d2 <- hill_derivative(x, order = 2, from = from, to = to, n = n,
                        domain = domain, allow_negative = allow_negative)
  abs(d2) / (1 + d1^2)^(3 / 2)
}

hill_radius <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                        allow_negative = FALSE) {
  curv <- hill_curvature(x, from = from, to = to, n = n, domain = domain,
                         allow_negative = allow_negative)
  out <- 1 / curv
  out[curv < 1e-8] <- NA_real_
  out
}

.area_vector <- function(y, baseline) {
  q <- attr(y, "q")
  .finite_grid(q, "Profile areas")
  out <- apply(y - baseline, 2L, function(value) .trapz(q, value))
  names(out) <- colnames(y)
  out
}

.arc_vector <- function(y) {
  q <- attr(y, "q")
  .finite_grid(q, "Arc lengths")
  out <- apply(y, 2L, function(value) sum(sqrt(diff(q)^2 + diff(value)^2)))
  names(out) <- colnames(y)
  out
}

.check_baseline <- function(baseline) {
  if (length(baseline) != 1L || !is.numeric(baseline) ||
      !is.finite(baseline)) {
    stop("'baseline' must be a single finite number.", call. = FALSE)
  }
  baseline
}

hill_area <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                      baseline = 0, allow_negative = FALSE) {
  baseline <- .check_baseline(baseline)
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  out <- .area_vector(y, baseline)
  data.frame(hill_area = unname(out), row.names = names(out))
}

hill_cum <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                     baseline = 0, allow_negative = FALSE) {
  baseline <- .check_baseline(baseline)
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  q <- attr(y, "q")
  .finite_grid(q, "Cumulative integrals")
  out <- apply(y - baseline, 2L, function(value) .cumtrapz(q, value))
  .as_profile(out, q, colnames(y))
}

hill_arc <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                     allow_negative = FALSE) {
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  out <- .arc_vector(y)
  data.frame(arc_length = unname(out), row.names = names(out))
}

hill_ranking <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                         baseline = 0, allow_negative = FALSE) {
  baseline <- .check_baseline(baseline)
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  aa <- .area_vector(y, baseline)
  ar <- .arc_vector(y)
  key <- .hill_matrix(x, c(0, 1, 2), which = 1L)
  data.frame(
    richness = key[1L, ],
    shannon = key[2L, ],
    simpson = key[3L, ],
    hill_area = unname(aa),
    area_rank = rank(-aa, ties.method = "average"),
    arc_length = unname(ar),
    arc_rank = rank(ar, ties.method = "average"),
    row.names = names(aa)
  )
}

hill_total_change <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                              baseline = 0, allow_negative = FALSE) {
  baseline <- .check_baseline(baseline)
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  aa <- .area_vector(y, baseline)
  k <- length(aa)
  if (k < 2L) {
    stop("At least two communities or time points are required.",
         call. = FALSE)
  }
  data.frame(
    from = names(aa)[-k],
    to = names(aa)[-1L],
    area_from = unname(aa[-k]),
    area_to = unname(aa[-1L]),
    absolute_change = unname(diff(aa)),
    relative_change = unname(diff(aa) / aa[-k]),
    stringsAsFactors = FALSE
  )
}
