x <- matrix(c(54, 20, 17,  6,  3,
              25, 25, 25, 25,  0,
              35, 35, 30,  0,  0),
            nrow = 3, byrow = TRUE,
            dimnames = list(c("A", "B", "C"), NULL))

reference <- function(p, q) {
  p <- p[p > 0]
  if (q == 1) exp(-sum(p * log(p))) else sum(p^q)^(1 / (1 - q))
}

test_that("profiles match known indices and the direct formula", {
  p <- x / rowSums(x)
  hp <- hill_profile(x, domain = c(0, 1, 2, Inf))
  expect_equal(unname(hp[1, ]), c(5, 4, 3))
  expect_equal(unname(hp[2, ]), unname(exp(-rowSums(ifelse(p > 0, p * log(p), 0)))))
  expect_equal(unname(hp[3, ]), unname(1 / rowSums(p^2)))
  expect_equal(unname(hp[4, ]), unname(1 / apply(p, 1, max)))

  q <- c(0.3, 0.999, 0.9995, 1, 1.0005, 1.001, 1.7, 9, 60)
  expected <- sapply(1:3, function(i) sapply(q, function(v) reference(p[i, ], v)))
  expect_equal(unname(hill_profile(x, domain = q)), expected,
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("profiles are non-increasing and flat for even communities", {
  hp <- hill_profile(x, from = 0, to = 30, n = 601)
  expect_true(all(diff(hp) <= 1e-10))
  expect_equal(unname(range(hp[, "B"])), c(4, 4))
  expect_equal(unname(hill_profile(x, domain = 1e5)[1, "A"]), 100 / 54,
               tolerance = 1e-4)
})

test_that("absolute and relative abundances give the same profile", {
  expect_equal(hill_profile(x), hill_profile(x / rowSums(x)))
})

test_that("single orders, vectors and data frames are accepted", {
  expect_equal(dim(hill_profile(x, domain = 2)), c(1L, 3L))
  expect_equal(unname(hill_profile(c(1, 1, 2), domain = 0)[1, 1]), 3)
  expect_equal(hill_profile(as.data.frame(x)), hill_profile(x),
               ignore_attr = TRUE)
  expect_equal(attr(hill_profile(x, from = 0, to = 10, n = 121), "q"),
               seq(0, 10, length.out = 121))
})

test_that("invalid input is rejected", {
  expect_error(hill_profile(x, domain = -1), "Negative orders")
  expect_error(hill_profile(-x))
  expect_error(hill_profile(rbind(x, 0)))
  expect_error(hill_profile(data.frame(a = "a", b = 1)))
  expect_error(hill_profile(x, from = 0, to = Inf))
  expect_error(hill_area(x, domain = c(2, 0, 1)), "strictly increasing")
  expect_error(hill_area(x, domain = c(0, Inf)), "finite")
  expect_error(hill_derivative(x, order = 3))
  expect_error(hill_total_change(x[1, , drop = FALSE]))
})

test_that("negative orders exceed richness when allowed", {
  neg <- hill_profile(x, domain = c(-2, -Inf), allow_negative = TRUE)
  expect_true(all(neg[1, ] >= c(5, 4, 3)))
  expect_equal(unname(neg[2, "A"]), 100 / 3)
})

test_that("analytical derivatives agree with finite differences", {
  q <- c(0.3, 0.9, 0.9995, 1, 1.0005, 1.3, 2, 6)
  h <- 1e-5
  d1 <- (hill_profile(x, domain = q + h) - hill_profile(x, domain = q - h)) /
    (2 * h)
  expect_equal(hill_first_derivative(x, domain = q), d1, tolerance = 1e-6,
               ignore_attr = TRUE)
  d2 <- (hill_first_derivative(x, domain = q + h) -
           hill_first_derivative(x, domain = q - h)) / (2 * h)
  expect_equal(hill_second_derivative(x, domain = q), d2, tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_true(all(hill_first_derivative(x) <= 1e-12))
})

test_that("curvature and radius are consistent", {
  d1 <- hill_first_derivative(x)
  d2 <- hill_second_derivative(x)
  k <- hill_curvature(x)
  expect_equal(k, abs(d2) / (1 + d1^2)^1.5)
  r <- hill_radius(x)
  expect_true(all(is.na(r[, "B"])))
  expect_equal(r[, "A"] * k[, "A"], rep(1, nrow(k)), ignore_attr = TRUE)
})

test_that("areas, cumulative integrals and arcs are coherent", {
  hc <- hill_cum(x, from = 0, to = 4, n = 401)
  ha <- hill_area(x, from = 0, to = 4, n = 401)
  expect_equal(unname(hc[1, ]), c(0, 0, 0))
  expect_equal(unname(hc[nrow(hc), ]), ha$hill_area)
  expect_equal(rownames(ha), c("A", "B", "C"))
  expect_equal(ha["B", 1], 16)
  expect_equal(hill_area(x, from = 0, to = 4, n = 401, baseline = 1)$hill_area,
               ha$hill_area - 4)
  p <- prop.table(x[1, ])
  exact <- integrate(function(q) sapply(q, function(v) reference(p, v)),
                     0, 4, rel.tol = 1e-10)$value
  expect_equal(ha["A", 1], exact, tolerance = 1e-5)
  arc <- hill_arc(x, from = 0, to = 4, n = 401)
  expect_equal(arc["B", 1], 4)
  expect_true(all(arc$arc_length >= 4))
})

test_that("ranking and change tables have the documented structure", {
  rk <- hill_ranking(x)
  expect_equal(rk$richness, c(5, 4, 3))
  expect_equal(rk$area_rank, rank(-rk$hill_area))
  expect_equal(rk$arc_rank[2], 1)
  tc <- hill_total_change(x)
  expect_equal(tc$from, c("A", "B"))
  expect_equal(tc$to, c("B", "C"))
  expect_equal(tc$relative_change, tc$absolute_change / tc$area_from)
  expect_equal(tc$area_to[1], tc$area_from[2])
})

test_that("surfaces use the expected grid and time", {
  s <- hill_surface(x, from = 0, to = 5, n = 11)
  expect_equal(dim(s$z), c(11L, 3L))
  expect_equal(s$time, 1:3)
  z <- x
  rownames(z) <- c(2001, 2005, 2010)
  expect_equal(hill_surface(z)$time, c(2001, 2005, 2010))
  expect_error(hill_surface(x, time = c(3, 2, 1)))
  expect_error(hill_surface(x[1, , drop = FALSE]))
  expect_error(hill_surface_compare(x, x[1:2, ], plot = "static"))
})

test_that("static plots run and return their data", {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })
  expect_equal(hill_profile_plot(x), hill_profile(x))
  expect_silent(hill_derivative_plot(x, order = 2))
  expect_silent(hill_radius_plot(x))
  s <- hill_surface_plot(x, plot = "static")
  expect_named(s, c("q", "time", "z"))
  cmp <- hill_surface_compare(x, x[3:1, ], plot = "static")
  expect_equal(cmp$difference$z, cmp$surface2$z - cmp$surface1$z)
  expect_silent(hill_surface_compare(x, x, plot = "static",
                                     mode = "difference"))
})

test_that("interactive plots return plotly objects", {
  skip_if_not_installed("plotly")
  expect_s3_class(hill_surface_plot(x), "plotly")
  expect_s3_class(hill_surface_compare(x, x[3:1, ]), "plotly")
  expect_s3_class(hill_surface_compare(x, x[3:1, ], mode = "difference"),
                  "plotly")
})

test_that("export_plots writes the requested file only", {
  expect_error(export_plots(x))
  f <- tempfile(fileext = ".pdf")
  expect_equal(export_plots(x, file = f), f)
  expect_true(file.exists(f))
  unlink(f)
})

test_that("names of BioFTF 1.x are deprecated and map to Hill tools", {
  expect_warning(out <- area(x), "deprecated")
  expect_equal(out, hill_area(x))
  expect_warning(out <- beta(x, n = 21), "deprecated")
  expect_equal(out, hill_profile(x, n = 21))
  expect_warning(out <- betasecond(x), "deprecated")
  expect_equal(out, hill_second_derivative(x))
  expect_warning(out <- ranking(x), "deprecated")
  expect_equal(out, hill_ranking(x))
})
