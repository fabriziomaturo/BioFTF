\name{datirel}
\alias{datirel}
\alias{summary_species}
\alias{summary_species_relative}
\alias{alltools}
\alias{export_plots}
\title{Data Summaries and Wrappers}
\description{
Relative abundances, species summaries, a wrapper that returns all the
functional tools, and a wrapper that saves the profile plots to a PDF file.
}
\usage{
datirel(x)
summary_species(x)
summary_species_relative(x)
alltools(x, from = 0, to = 2, n = 101)
export_plots(x, file, from = 0, to = 2, n = 101)
}
\arguments{
\item{x}{A numeric matrix or data frame of non-negative abundances. Rows are
communities, sites, or time points; columns are species.}
\item{from, to}{Finite lower and upper endpoints of the grid of orders.}
\item{n}{Number of equally spaced orders.}
\item{file}{Path of the PDF file to create.}
}
\value{
\code{datirel} returns the matrix of relative abundances.
\code{summary_species} returns a data frame with richness, total abundance,
and smallest and largest abundance of each community, zeros included;
\code{summary_species_relative} returns the same summaries for relative
abundances. \code{alltools} returns a list with the relative abundances, the
species summary, and the output of \code{\link{hill_profile}},
\code{\link{hill_first_derivative}}, \code{\link{hill_second_derivative}},
\code{\link{hill_curvature}}, \code{\link{hill_radius}},
\code{\link{hill_area}}, \code{\link{hill_cum}}, \code{\link{hill_arc}}, and
\code{\link{hill_ranking}}. \code{export_plots} writes one page for each of
profile, first derivative, second derivative, curvature, radius of
curvature, and cumulative integral function, and invisibly returns
\code{file}.
}
\examples{
x <- matrix(c(3, 5, 5, 2, 1,
              2, 6, 8, 5, 1,
              7, 3, 2, 1, 0),
            nrow = 3, byrow = TRUE)
datirel(x)
summary_species(x)
summary_species_relative(x)

res <- alltools(x, from = 0, to = 5, n = 251)
names(res)
res$ranking

out <- file.path(tempdir(), "BioFTF-plots.pdf")
export_plots(x, file = out)
unlink(out)
}
\keyword{multivariate}
