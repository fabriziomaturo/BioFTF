\name{hill_surface}
\alias{hill_surface}
\alias{hill_surface_plot}
\alias{hill_surface_compare}
\title{Hill-Number Biodiversity Surfaces}
\description{
Compute, plot, and compare surfaces obtained by placing the Hill-number
profiles of ordered rows side by side. Plots are either rotatable or static.
}
\usage{
hill_surface(x, from = 0, to = 2, n = 101, domain = NULL,
             time = NULL, allow_negative = FALSE)
hill_surface_plot(x, from = 0, to = 2, n = 101, domain = NULL,
                  time = NULL, plot = c("interactive", "static"),
                  allow_negative = FALSE,
                  xlab = "Order q", ylab = "Time",
                  zlab = "Hill number",
                  main = "Hill-number surface",
                  theta = 40, phi = 25, ...)
hill_surface_compare(x, y, from = 0, to = 2, n = 101,
                     domain = NULL, time = NULL,
                     labels = c("Surface 1", "Surface 2"),
                     plot = c("interactive", "static"),
                     mode = c("overlay", "difference"),
                     allow_negative = FALSE,
                     main = "Comparison of Hill-number surfaces",
                     theta = 40, phi = 25, ...)
}
\arguments{
\item{x, y}{Numeric abundance matrices or data frames with at least two rows.
Rows are ordered time points or sites along a gradient; columns are species.
\code{x} and \code{y} must have the same number of rows; the species may
differ.}
\item{from, to}{Finite lower and upper endpoints of the grid of orders.}
\item{n}{Number of equally spaced orders used when \code{domain = NULL}.}
\item{domain}{Optional vector of finite, strictly increasing orders.}
\item{time}{Optional numeric vector of strictly increasing values, one for
each row. If \code{NULL}, row names are used when they can be read as
increasing numbers, for example years; otherwise rows are numbered from 1.}
\item{plot}{\code{"interactive"} (default) returns a \pkg{plotly} surface
that can be rotated and zoomed in the viewer or in a browser.
\code{"static"} draws a perspective plot on the current graphics device.}
\item{allow_negative}{Logical. If \code{FALSE} (default), negative orders
raise an error. See \code{\link{hill_profile}}.}
\item{xlab, ylab, zlab, main}{Axis labels and title.}
\item{theta, phi}{Viewing angles of the static plot, as in
\code{\link[graphics]{persp}}.}
\item{labels}{Names of the two surfaces.}
\item{mode}{\code{"overlay"} shows the two surfaces together (interactive) or
side by side on a common vertical scale (static). \code{"difference"} shows
the pointwise difference \code{y} minus \code{x}.}
\item{...}{Further arguments passed to \code{\link[plotly]{plot_ly}} for
interactive plots or to \code{\link[graphics]{persp}} for static plots.}
}
\details{
The surface has the order \eqn{q} on one horizontal axis, time on the other,
and the Hill number on the vertical axis. Profiles are observed only at the
rows of the data: between two consecutive rows the surface is a linear
interpolation drawn for display, not an estimate of diversity at
intermediate times.

Interactive plots require the suggested package \pkg{plotly}. In the static
plot the colour of each facet follows the height of the surface, and lines
mark the observed profiles and a subset of the orders.
}
\value{
\code{hill_surface} returns a list with components \code{q} (orders),
\code{time}, and \code{z}, a matrix of Hill numbers with orders in rows and
time points in columns. With \code{plot = "interactive"},
\code{hill_surface_plot} and \code{hill_surface_compare} return a
\pkg{plotly} object. With \code{plot = "static"}, \code{hill_surface_plot}
invisibly returns the surface, and \code{hill_surface_compare} invisibly
returns a list with the two surfaces and their difference.
}
\references{
Maturo, F. and Di Battista, T. (2018). A functional approach to Hill's numbers
for assessing changes in species variety of ecological communities over time.
\emph{Ecological Indicators}, 84, 70--81.
\doi{10.1016/j.ecolind.2017.08.016}
}
\seealso{
\code{\link{hill_profile}}, \code{\link{hill_total_change}}.
}
\examples{
years <- 2018:2027
x <- matrix(c(40, 35, 15,  5,  5,  0,
              38, 34, 16,  6,  4,  2,
              35, 30, 18,  8,  5,  4,
              30, 28, 20, 10,  7,  5,
              25, 26, 21, 12,  9,  7,
              20, 24, 22, 14, 11,  9,
              16, 21, 23, 16, 13, 11,
              12, 18, 24, 18, 15, 13,
               9, 15, 24, 20, 17, 15,
               6, 12, 23, 22, 19, 18),
            nrow = length(years), byrow = TRUE)
rownames(x) <- years

## A second set of surveys on the same years.
y <- x
y[, 1] <- x[, 1] + 15
y[, 6] <- x[, 6] + 4

## Surface between q = 0 and q = 10.
s <- hill_surface(x, from = 0, to = 10, n = 121)
str(s)

## Static figures.
hill_surface_plot(x, from = 0, to = 10, n = 121, plot = "static")
hill_surface_plot(x, from = 0, to = 10, n = 121, plot = "static",
                  theta = 130, phi = 20,
                  main = "Hill-number surface, 2018-2027")
hill_surface_compare(x, y, from = 0, to = 10, n = 121,
                     labels = c("Area 1", "Area 2"), plot = "static")
hill_surface_compare(x, y, from = 0, to = 10, n = 121,
                     labels = c("Area 1", "Area 2"), plot = "static",
                     mode = "difference")

## Rotatable figures.
if (interactive() && requireNamespace("plotly", quietly = TRUE)) {
  hill_surface_plot(x, from = 0, to = 10, n = 121)
  hill_surface_compare(x, y, from = 0, to = 10, n = 121,
                       labels = c("Area 1", "Area 2"))
  hill_surface_compare(x, y, from = 0, to = 10, n = 121,
                       labels = c("Area 1", "Area 2"),
                       mode = "difference")
}
}
\keyword{multivariate}
\keyword{hplot}
