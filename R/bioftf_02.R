## ---------------------------------------------------------------------------
## 2. Hill number and its analytical derivatives
## ---------------------------------------------------------------------------

.hill_core <- function(p, q) {
  p <- p[p > 0]
  lp <- log(p)
  e <- q - 1
  if (abs(e) < 1e-3) {
    m1 <- sum(p * lp)
    d <- lp - m1
    mu2 <- sum(p * d^2)
    mu3 <- sum(p * d^3)
    mu4 <- sum(p * d^4)
    mu5 <- sum(p * d^5)
    mu6 <- sum(p * d^6)
    k4 <- mu4 - 3 * mu2^2
    k5 <- mu5 - 10 * mu3 * mu2
    k6 <- mu6 - 15 * mu4 * mu2 - 10 * mu3^2 + 30 * mu2^3
    g <- -(m1 + mu2 * e / 2 + mu3 * e^2 / 6 + k4 * e^3 / 24 +
             k5 * e^4 / 120 + k6 * e^5 / 720)
    g1 <- -(mu2 / 2 + mu3 * e / 3 + k4 * e^2 / 8 + k5 * e^3 / 30 +
              k6 * e^4 / 144)
    g2 <- -(mu3 / 3 + k4 * e / 4 + k5 * e^2 / 10 + k6 * e^3 / 36)
  } else {
    a <- q * lp
    amax <- max(a)
    w <- exp(a - amax)
    sw <- sum(w)
    w <- w / sw
    L <- if (abs(e) < 0.5) log1p(sum(p * expm1(e * lp))) else amax + log(sw)
    m <- sum(w * lp)
    v <- sum(w * (lp - m)^2)
    u <- 1 - q
    g <- L / u
    g1 <- m / u + L / u^2
    g2 <- v / u + 2 * m / u^2 + 2 * L / u^3
  }
  D <- exp(g)
  c(D, D * g1, D * (g1^2 + g2))
}

.hill_number <- function(p, q) {
  p <- p[p > 0]
  if (is.infinite(q)) {
    return(if (q > 0) 1 / max(p) else 1 / min(p))
  }
  if (q == 0) {
    return(length(p))
  }
  .hill_core(p, q)[1L]
}

.hill_matrix <- function(x, q, which = 1L) {
  p <- .relative_abundance(x)
  values <- vapply(seq_len(nrow(p)), function(i) {
    vapply(q, function(value) {
      if (which == 1L) .hill_number(p[i, ], value)
      else .hill_core(p[i, ], value)[which]
    }, numeric(1L))
  }, numeric(length(q)))
  .as_profile(values, q, rownames(p))
}
