## ---------------------------------------------------------------------------
## 1. Input checks and grids
## ---------------------------------------------------------------------------

.check_abundance_matrix <- function(x) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1L)))) {
      stop("'x' must contain numeric columns only.", call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (is.numeric(x) && is.null(dim(x))) {
    x <- matrix(x, nrow = 1L, dimnames = list(NULL, names(x)))
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'x' must be a numeric matrix or data frame.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (anyNA(x)) {
    stop("'x' must not contain missing values.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("'x' must contain only finite values.", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("'x' must contain non-negative abundances.", call. = FALSE)
  }
  if (any(rowSums(x) <= 0)) {
    stop("Each community must have a positive total abundance.", call. = FALSE)
  }
  if (is.null(rownames(x))) {
    rownames(x) <- paste("community", seq_len(nrow(x)), sep = ".")
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste("species", seq_len(ncol(x)), sep = ".")
  }
  x
}

.relative_abundance <- function(x) {
  x <- .check_abundance_matrix(x)
  x / rowSums(x)
}

.hill_domain <- function(domain = NULL, from = 0, to = 2, n = 101,
                         allow_negative = FALSE) {
  if (is.null(domain)) {
    if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 3) {
      stop("'n' must be a single number greater than or equal to 3.",
           call. = FALSE)
    }
    if (length(from) != 1L || length(to) != 1L ||
        !is.finite(from) || !is.finite(to)) {
      stop("'from' and 'to' must be finite numbers. Use 'domain' to ",
           "evaluate the profile at Inf.", call. = FALSE)
    }
    if (to <= from) {
      stop("'to' must be greater than 'from'.", call. = FALSE)
    }
    domain <- seq(from, to, length.out = as.integer(round(n)))
  }
  if (!is.numeric(domain) || length(domain) < 1L || anyNA(domain)) {
    stop("'domain' must be a numeric vector without missing values.",
         call. = FALSE)
  }
  domain <- as.numeric(domain)
  if (any(domain < 0) && !allow_negative) {
    stop("Negative orders are rejected by default: for q < 0 the Hill ",
         "number exceeds species richness and is driven by the rarest ",
         "species. Set 'allow_negative = TRUE' to compute it anyway.",
         call. = FALSE)
  }
  domain
}

## Grid required by tools that integrate or differentiate along q.
.finite_grid <- function(q, what) {
  if (any(!is.finite(q))) {
    stop(what, " require finite orders.", call. = FALSE)
  }
  if (length(q) < 2L || any(diff(q) <= 0)) {
    stop(what, " require at least two strictly increasing orders.",
         call. = FALSE)
  }
  invisible(q)
}

.q_labels <- function(q) format(q, digits = 10, trim = TRUE)

.as_profile <- function(values, q, communities) {
  out <- matrix(values, nrow = length(q), ncol = length(communities),
                dimnames = list(.q_labels(q), communities))
  attr(out, "q") <- q
  out
}

.trapz <- function(x, y) {
  sum(diff(x) * (y[-length(y)] + y[-1L]) / 2)
}

.cumtrapz <- function(x, y) {
  c(0, cumsum(diff(x) * (y[-length(y)] + y[-1L]) / 2))
}

