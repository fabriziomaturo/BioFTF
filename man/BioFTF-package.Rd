\name{BioFTF-package}
\alias{BioFTF-package}
\alias{BioFTF}
\docType{package}
\title{Biodiversity Assessment Using Hill-Number Functional Tools}
\description{
Hill-number diversity profiles are treated as functions of the order
\eqn{q} and summarised through functional tools: derivatives, curvature,
radius of curvature, arc length, cumulative integral functions, areas,
rankings, and biodiversity surfaces over time.
}
\details{
Rows of the input matrix are communities, sites, or time points. Columns are
species. Absolute and relative non-negative abundances are both accepted,
because every row is rescaled to sum to one.

The main functions are \code{\link{hill_profile}} for the profiles,
\code{\link{hill_derivative}}, \code{\link{hill_curvature}},
\code{\link{hill_area}}, \code{\link{hill_cum}}, \code{\link{hill_arc}} and
\code{\link{hill_ranking}} for the functional tools, and
\code{\link{hill_surface_plot}} and \code{\link{hill_surface_compare}} for
static or rotatable three-dimensional surfaces.

Versions 1.x of the package computed the same functional tools on the
Patil-Taillie beta profile. The relation between the two profiles is given in
\code{\link{hill_profile}}. The function names of versions 1.x are still
available but deprecated; see \code{\link{BioFTF-deprecated}}.
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

Hill, M. O. (1973). Diversity and evenness: a unifying notation and its
consequences. \emph{Ecology}, 54(2), 427--432. \doi{10.2307/1934352}

Patil, G. P. and Taillie, C. (1982). Diversity as a concept and its
measurement. \emph{Journal of the American Statistical Association}, 77(379),
548--561. \doi{10.1080/01621459.1982.10477845}
}
\author{
Fabrizio Maturo, Francesca Fortuna, Tonio Di Battista.
}
\keyword{package}
