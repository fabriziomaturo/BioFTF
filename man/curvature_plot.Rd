\name{curvature_plot}
\alias{curvature_plot}
\title{Plot the curvature of the beta profiles for each community.}
\usage{
curvature_plot(x)
}
\arguments{\item{x}{A data matrix with the abundance of the species (the rows are the communities and the columns indicate the species living in a community).}}

\description{
This function provides a plot of the curvature functions of the diversity profiles for the considered communities. The domain "beta" represents the sensibility to changes in the presence of rare species. The graph of the curvature reflects the community composition; in particular, if we consider communities with different numbers of species, the curvature is a good indicator of the lack of balance among species. Great levels of curvature highlights high imbalance among different species and vice-versa.
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
curvature_plot(x)
}
