\name{area}
\alias{area}
\title{Compute the area under the beta profile for each community.}
\usage{
area(x)
}
\arguments{\item{x}{A data matrix with the abundance of the species (the rows are the communities and the columns indicate the species living in a community).}}

\description{
This function provides the area under the beta profile for each community. In an ecological framework, the area under the profiles can be used to assess biodiversity because it can rank communities, even if they have different number of species. Indeed, in a case of maximum dominance, the area is extremely low; while in a case of evenness, it increase. The main advantage of using area is that the ordering among communities can be investigated without the analysis of a graph. It provides a scalar measure of biodiversity preserving its multidimensional aspect. This is not possible with the classical indices.
}
\examples{
x=matrix(c(0.3,0.5,0.1,0.05,0.05,0.25,0.25,0.25,0.25,0,0.35,0.3,0.35,0,0),3,5)
area(x)

}
