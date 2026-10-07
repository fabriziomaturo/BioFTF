\name{BioFTF-deprecated}
\alias{BioFTF-deprecated}
\alias{beta}
\alias{beta_plot}
\alias{betaprime}
\alias{betaprime_plot}
\alias{betasecond}
\alias{betasecond_plot}
\alias{curvature}
\alias{curvature_plot}
\alias{radius}
\alias{radius_plot}
\alias{arc}
\alias{area}
\alias{ranking}
\title{Deprecated Functions}
\description{
Function names of BioFTF 1.x, kept so that existing code still runs. Each
function issues a warning and calls its replacement.
}
\usage{
beta(x, n = 101, ...)
beta_plot(x, n = 101, ...)
betaprime(x, n = 101, ...)
betaprime_plot(x, n = 101, ...)
betasecond(x, n = 101, ...)
betasecond_plot(x, n = 101, ...)
curvature(x, n = 101, ...)
curvature_plot(x, n = 101, ...)
radius(x, n = 101, ...)
radius_plot(x, n = 101, ...)
arc(x, n = 101, ...)
area(x, n = 101, ...)
ranking(x, n = 101, ...)
}
\arguments{
\item{x}{A numeric matrix or data frame of non-negative abundances.}
\item{n}{Number of equally spaced orders.}
\item{...}{Further arguments passed to the replacement function.}
}
\details{
In BioFTF 1.x these functions were computed on the Patil-Taillie beta
profile over \eqn{-1 \le \beta \le 1}{-1 <= beta <= 1}. They now return the
corresponding tool for the Hill-number profile over \eqn{0 \le q \le 2}, so
their values differ from those of BioFTF 1.x. The relation between the two
profiles is described in \code{\link{hill_profile}}.

\tabular{ll}{
Deprecated \tab Replacement \cr
\code{beta}, \code{beta_plot} \tab \code{\link{hill_profile}}, \code{\link{hill_profile_plot}} \cr
\code{betaprime}, \code{betaprime_plot} \tab \code{\link{hill_first_derivative}}, \code{\link{hill_derivative_plot}} \cr
\code{betasecond}, \code{betasecond_plot} \tab \code{\link{hill_second_derivative}}, \code{\link{hill_derivative_plot}} \cr
\code{curvature}, \code{curvature_plot} \tab \code{\link{hill_curvature}}, \code{\link{hill_curvature_plot}} \cr
\code{radius}, \code{radius_plot} \tab \code{\link{hill_radius}}, \code{\link{hill_radius_plot}} \cr
\code{arc} \tab \code{\link{hill_arc}} \cr
\code{area} \tab \code{\link{hill_area}} \cr
\code{ranking} \tab \code{\link{hill_ranking}} \cr
}

\code{beta} masks \code{\link[base:Special]{beta}} of package \pkg{base}
while BioFTF is attached; the Beta function remains available as
\code{base::beta}.
}
\value{
The value of the replacement function.
}
\seealso{
\code{\link{hill_profile}}, \code{\link{hill_derivative}}.
}
\keyword{internal}
