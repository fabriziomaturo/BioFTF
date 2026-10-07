\name{hill_derivative}
\alias{hill_derivative}
\alias{hill_first_derivative}
\alias{hill_second_derivative}
\alias{hill_curvature}
\alias{hill_radius}
\alias{hill_area}
\alias{hill_cum}
\alias{hill_arc}
\alias{hill_ranking}
\alias{hill_total_change}
\alias{hill_derivative_plot}
\alias{hill_curvature_plot}
\alias{hill_radius_plot}
\alias{hill_cum_plot}
\title{Functional Tools Derived from Hill-Number Profiles}
\description{
Derivatives, curvature, radius of curvature, cumulative integral functions,
areas, arc lengths, rankings, and changes between consecutive rows, all
computed from Hill-number profiles.
}
\usage{
hill_derivative(x, order = 1, from = 0, to = 2, n = 101,
                domain = NULL, allow_negative = FALSE)
hill_first_derivative(x, from = 0, to = 2, n = 101,
                      domain = NULL, allow_negative = FALSE)
hill_second_derivative(x, from = 0, to = 2, n = 101,
                       domain = NULL, allow_negative = FALSE)
hill_curvature(x, from = 0, to = 2, n = 101,
               domain = NULL, allow_negative = FALSE)
hill_radius(x, from = 0, to = 2, n = 101,
            domain = NULL, allow_negative = FALSE)
hill_area(x, from = 0, to = 2, n = 101, domain = NULL,
          baseline = 0, allow_negative = FALSE)
hill_cum(x, from = 0, to = 2, n = 101, domain = NULL,
         baseline = 0, allow_negative = FALSE)
hill_arc(x, from = 0, to = 2, n = 101,
         domain = NULL, allow_negative = FALSE)
hill_ranking(x, from = 0, to = 2, n = 101, domain = NULL,
             baseline = 0, allow_negative = FALSE)
hill_total_change(x, from = 0, to = 2, n = 101, domain = NULL,
                  baseline = 0, allow_negative = FALSE)
hill_derivative_plot(x, order = 1, from = 0, to = 2, n = 101,
                     domain = NULL, allow_negative = FALSE,
                     legend = "bottomright", ...)
hill_curvature_plot(x, from = 0, to = 2, n = 101,
                    domain = NULL, allow_negative = FALSE,
                    legend = "topright", ...)
hill_radius_plot(x, from = 0, to = 2, n = 101,
                 domain = NULL, allow_negative = FALSE,
                 legend = "topleft", ...)
hill_cum_plot(x, from = 0, to = 2, n = 101, domain = NULL,
              baseline = 0, allow_negative = FALSE,
              legend = "topleft", ...)
}
\arguments{
\item{x}{A numeric matrix or data frame of non-negative abundances. Rows are
communities, sites, or time points; columns are species.}
\item{order}{Derivative order, 1 or 2.}
\item{from, to}{Finite lower and upper endpoints of the grid of orders.}
\item{n}{Number of equally spaced orders used when \code{domain = NULL}.}
\item{domain}{Optional numeric vector of finite orders. Areas, cumulative
integrals, and arc lengths require at least two strictly increasing values.}
\item{baseline}{Constant subtracted from the profile before integration.
With \code{baseline = 1}, the smallest value a Hill number can take, the area
of a community with a single species is zero.}
\item{allow_negative}{Logical. If \code{FALSE} (default), negative orders
raise an error. See \code{\link{hill_profile}}.}
\item{legend}{Position of the legend, as accepted by
\code{\link[graphics]{legend}}, or \code{NULL} for no legend.}
\item{...}{Further arguments passed to \code{\link[graphics]{matplot}}.
\code{hill_radius_plot} uses a logarithmic vertical axis.}
}
\details{
\emph{Derivatives.} First and second derivatives are analytical. Writing
\eqn{L(q) = \log \sum_i p_i^q}{L(q) = log(sum_i p_i^q)} and
\eqn{g(q) = L(q)/(1 - q)}, the Hill number is \eqn{D = \exp(g)}{D = exp(g)},
so that \eqn{D' = D g'} and \eqn{D'' = D (g'^2 + g'')}. The derivatives of
\eqn{L} are the mean and the variance of \eqn{\log p_i}{log(p_i)} under the
weights \eqn{p_i^q / \sum_j p_j^q}{p_i^q / sum_j p_j^q}. Near \eqn{q = 1},
where the closed form is numerically unstable, a Taylor expansion is used.
The derivatives therefore do not depend on the spacing of the grid.

\emph{Curvature and radius.} The curvature is
\eqn{|D''| / (1 + D'^2)^{3/2}} and the radius of curvature is its
reciprocal. Where the curvature is numerically zero (below
\eqn{10^{-8}}), as for a community whose species are equally abundant or
far along the flat tail of a profile, the radius is returned as \code{NA}.

\emph{Area and cumulative integral.} \code{hill_cum} returns
\eqn{\int_{q_1}^{q} (D(u) - b)\,du}{the integral of D(u) - b from q_1 to q},
with \eqn{b} the baseline, at each order of the grid, and \code{hill_area}
its value at the last order. Both use the trapezoidal rule, whose error
decreases with the square of the grid spacing. With the default
\code{n = 101} the relative error of the area is of order \eqn{10^{-4}} on
\eqn{[0, 2]} and grows with wider domains, so \code{n} should be increased
together with \code{to}.

\emph{Arc length.} \code{hill_arc} returns the length of the polygonal line
through the profile. It is bounded below by \code{to - from}, attained by a
flat profile. The arc length depends on the units of both axes.

\emph{Ranking.} \code{hill_ranking} reports the Hill numbers of order 0, 1,
and 2 together with area and arc length. \code{area_rank} is 1 for the
largest area. \code{arc_rank} is 1 for the shortest arc, that is, for the
flattest profile; it orders communities by how little the profile changes
along the domain and is not a ranking by diversity. A larger area does not
imply that one profile lies above the other at every order: when two
profiles cross, the communities are not ordered in an intrinsic sense and the
area ranking depends on the chosen domain. As \code{to} increases, the area
is progressively dominated by \eqn{D(\infty) = 1/\max_i p_i}{D(Inf) = 1/max_i p_i}.

\emph{Change.} \code{hill_total_change} compares the area of each row with
the area of the previous one, which is meaningful when rows are ordered, for
instance in time.
}
\value{
\code{hill_derivative}, \code{hill_first_derivative},
\code{hill_second_derivative}, \code{hill_curvature}, \code{hill_radius}, and
\code{hill_cum} return matrices with orders in rows and communities in
columns; the exact orders are stored in the attribute \code{"q"}.
\code{hill_area} and \code{hill_arc} return one-column data frames with one
row per community. \code{hill_ranking} returns a data frame with columns
\code{richness}, \code{shannon}, \code{simpson}, \code{hill_area},
\code{area_rank}, \code{arc_length}, and \code{arc_rank}.
\code{hill_total_change} returns a data frame with columns \code{from},
\code{to}, \code{area_from}, \code{area_to}, \code{absolute_change}, and
\code{relative_change}. Plot functions return the plotted matrix invisibly.
}
\references{
Di Battista, T., Fortuna, F. and Maturo, F. (2017). BioFTF: An R package for
biodiversity assessment with the functional data analysis approach.
\emph{Ecological Indicators}, 73, 726--732.
\doi{10.1016/j.ecolind.2016.10.032}

Maturo, F. and Di Battista, T. (2018). A functional approach to Hill's numbers
for assessing changes in species variety of ecological communities over time.
\emph{Ecological Indicators}, 84, 70--81.
\doi{10.1016/j.ecolind.2017.08.016}

Maturo, F. (2018). Unsupervised classification of ecological communities ranked
according to their biodiversity patterns via a functional principal component
decomposition of Hill's numbers integral functions.
\emph{Ecological Indicators}, 90, 305--315.
\doi{10.1016/j.ecolind.2018.03.013}
}
\seealso{
\code{\link{hill_profile}}, \code{\link{hill_surface_plot}}.
}
\examples{
x <- matrix(c(54, 20, 17,  6,  3,
              25, 25, 25, 25,  0,
              35, 35, 30,  0,  0),
            nrow = 3, byrow = TRUE)
rownames(x) <- c("A", "B", "C")

## Derivatives at selected orders.
hill_first_derivative(x, domain = c(0, 0.5, 1, 2))
hill_second_derivative(x, domain = c(0, 0.5, 1, 2))
head(hill_curvature(x))

## Area, cumulative integral function, arc length, and ranking on [0, 2].
hill_area(x)
tail(hill_cum(x))
hill_arc(x)
hill_ranking(x)

## Area above the minimum value of a Hill number on a wider domain.
hill_area(x, from = 0, to = 10, n = 1001, baseline = 1)

## Plots.
hill_derivative_plot(x, order = 1)
hill_curvature_plot(x)
hill_cum_plot(x, from = 0, to = 10, n = 501)

## Change in area between consecutive surveys.
surveys <- matrix(c(40, 30, 15, 10,  5,
                    35, 30, 18, 12,  5,
                    30, 28, 20, 15,  7,
                    25, 25, 22, 18, 10,
                    20, 24, 24, 20, 12,
                    15, 22, 25, 22, 16),
                  nrow = 6, byrow = TRUE)
rownames(surveys) <- 2018:2023
hill_total_change(surveys)
}
\keyword{multivariate}
\keyword{hplot}
