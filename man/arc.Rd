\name{arc}
\alias{arc}
\title{Compute the arc length of beta profile for each community.}
\usage{
arc(x)
}
\arguments{\item{x}{A data matrix with the abundance of the species (the rows are the communities and the columns indicate the species living in a community).}}

\description{
This function provides the arc length of beta profile for each community. In an ecological framework, the arc length can be used to assess biodiversity because it can rank communities, given the number of species. Indeed, in a case of maximum dominance, the profile length is extremely high; while in a case of evenness, it decreases. The main advantage of using profile length is that the ordering among communities can be investigated without the analysis of a graph. It provides a scalar measure of diversity preserving its multidimensional aspect. This is not possible with the classical indices.
}

\author{
Fabrizio Maturo, Francesca Fortuna, Tonio Di Battista
}

\source{
ENVIRONMENTAL MONITORING THROUGH FUNCTIONAL BIODIVERSITY TOOLS - T. Di Battista, F. Fortuna, F. Maturo - Ecological Indicators (2016) 60, 237-247 - online (2015). DOI: 10.1016/j.ecolind.2015.05.056.

PARAMETRIC FUNCTIONAL ANALYSIS OF VARIANCE FOR FISH BIODIVERSITY ASSESSMENT - T. Di Battista, F. Fortuna, F. Maturo - Journal of Environmental Informatics. Special issue of Marine and Freshwater Environmental Engineering and Management, Journal of Environmental Informatics. In press.

}

\examples{
x=matrix(c(0.3,0.5,0.1,0.05,0.05,0.25,0.25,0.25,0.25,0,0.35,0.3,0.35,0,0),3,5)
arc(x)

}
