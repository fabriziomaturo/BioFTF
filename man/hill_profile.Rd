\name{hill_profile}
\alias{hill_profile}
\alias{hill}
\alias{hill_profile_plot}
\alias{hill_plot}
\title{Hill-Number Diversity Profiles}
\description{
Compute and plot Hill-number diversity profiles for ecological communities,
sites, or time points.
}
\usage{
hill_profile(x, from = 0, to = 2, n = 101, domain = NULL,
             allow_negative = FALSE)
hill(x, from = 0, to = 2, n = 101, domain = NULL,
     allow_negative = FALSE)
hill_profile_plot(x, from = 0, to = 2, n = 101, domain = NULL,
                  allow_negative = FALSE, legend = "topright", ...)
hill_plot(x, from = 0, to = 2, n = 101, domain = NULL,
          allow_negative = FALSE, legend = "topright", ...)
}
\arguments{
\item{x}{A numeric matrix or data frame of non-negative abundances. Rows are
communities, sites, or time points; columns are species. A numeric vector is
treated as a single community.}
\item{from, to}{Finite lower and upper endpoints of the grid of orders. The
default is \eqn{0 \le q \le 2}. Larger upper endpoints, such as 10, 20, or
100, extend the profile towards the dominance end.}
\item{n}{Number of equally spaced orders between \code{from} and \code{to},
used when \code{domain = NULL}.}
\item{domain}{Optional numeric vector of orders at which the profile is
evaluated. When supplied, \code{from}, \code{to} and \code{n} are ignored. It
may contain a single order and may contain \code{Inf}.}
\item{allow_negative}{Logical. If \code{FALSE} (default), negative orders
raise an error. See Details.}
\item{legend}{Position of the legend, as accepted by
\code{\link[graphics]{legend}}, or \code{NULL} for no legend.}
\item{...}{Further arguments passed to \code{\link[graphics]{matplot}}.}
}
\details{
For relative abundances \eqn{p_1, \dots, p_S}, the Hill number of order
\eqn{q} is
\deqn{D(q) = \left(\sum_{i: p_i > 0} p_i^q\right)^{1/(1-q)}, \qquad q \ne 1,}{D(q) = (sum_i p_i^q)^(1/(1-q)), q != 1,}
with \eqn{D(1) = \exp(-\sum_i p_i \log p_i)}{D(1) = exp(-sum_i p_i log p_i)}
defined by continuity. \eqn{D(0)} is species richness, \eqn{D(1)} is the
exponential of Shannon entropy, \eqn{D(2)} is the inverse Simpson
concentration, and \eqn{D(\infty) = 1/\max_i p_i}{D(Inf) = 1/max_i p_i} is the
reciprocal of the Berger-Parker dominance. \eqn{D(q)} is expressed in
effective numbers of species and is non-increasing in \eqn{q}; it is constant
only when all the species present are equally abundant.

\emph{Range of the order.} Most of the variation of the profile lies between
\eqn{q = 0} and \eqn{q = 2}, which is the default. Increasing \code{to} adds
the part of the curve that converges to \eqn{D(\infty)}{D(Inf)}; the limit
itself is obtained with \code{domain = Inf}.

\emph{Negative orders.} For \eqn{q < 0} the Hill number is larger than
species richness, so it can no longer be read as an effective number of
species out of the \eqn{S} observed, and it is driven by the smallest
positive abundances: \eqn{D(-\infty) = 1/\min_i p_i}{D(-Inf) = 1/min_i p_i}.
These are the abundances most affected by sampling error, and the value
changes abruptly when a rare species goes undetected. Negative orders are
therefore rejected unless \code{allow_negative = TRUE}.

\emph{Relation with the Patil-Taillie profile.} Versions 1.x of the package
used the beta profile of Patil and Taillie (1982),
\deqn{\Delta(\beta) = \frac{1 - \sum_i p_i^{\beta + 1}}{\beta}, \qquad -1 \le \beta \le 1,}{Delta(beta) = (1 - sum_i p_i^(beta + 1)) / beta, -1 <= beta <= 1,}
which equals \eqn{S - 1} at \eqn{\beta = -1}{beta = -1}, Shannon entropy as
\eqn{\beta \to 0}{beta -> 0}, and the Gini-Simpson index at
\eqn{\beta = 1}{beta = 1}. With \eqn{q = \beta + 1}{q = beta + 1} the two
profiles are linked by
\deqn{\Delta(\beta) = \frac{1 - D(q)^{-\beta}}{\beta},}{Delta(beta) = (1 - D(q)^(-beta)) / beta,}
so the interval \eqn{[-1, 1]} for \eqn{\beta}{beta} corresponds to
\eqn{[0, 2]} for \eqn{q}. The two profiles induce the same ordering of
communities at each order, but they are measured on different scales: the
beta profile is an entropy, whereas the Hill number is its numbers
equivalent. Derivatives, curvature, arc length, and areas are not invariant
to this transformation, so functional tools computed on the two profiles
differ and may rank communities differently.
}
\value{
\code{hill_profile} and \code{hill} return a matrix with orders in rows and
communities in columns. The exact orders are stored in the attribute
\code{"q"}. Plot functions return the profile invisibly.
}
\references{
Hill, M. O. (1973). Diversity and evenness: a unifying notation and its
consequences. \emph{Ecology}, 54(2), 427--432. \doi{10.2307/1934352}

Patil, G. P. and Taillie, C. (1982). Diversity as a concept and its
measurement. \emph{Journal of the American Statistical Association}, 77(379),
548--561. \doi{10.1080/01621459.1982.10477845}

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
\code{\link{hill_derivative}}, \code{\link{hill_area}},
\code{\link{hill_surface_plot}}.
}
\examples{
x <- matrix(c(54, 20, 17,  6,  3,
              25, 25, 25, 25,  0,
              35, 35, 30,  0,  0,
              70, 10, 10,  5,  5),
            nrow = 4, byrow = TRUE)
rownames(x) <- c("A", "B", "C", "D")

## Default grid: 101 orders between 0 and 2.
p <- hill_profile(x)
head(p)

## Selected orders, including the limit for q tending to infinity.
hill_profile(x, domain = c(0, 1, 2, 5, 10, 100, Inf))

## A wider domain.
hill_profile_plot(x, from = 0, to = 10, n = 201)

## Negative orders exceed species richness.
hill_profile(x, domain = c(-1, 0), allow_negative = TRUE)
}
\keyword{multivariate}
\keyword{hplot}
