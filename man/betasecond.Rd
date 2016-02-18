\name{betasecond}
\alias{betasecond}
\title{Compute the values of the second derivative of beta profile for each community in every part of the beta domain.}
\usage{
betasecond(x)
}
\arguments{\item{x}{A data matrix with the abundance of the species (the rows are the communities and the columns indicate the species living in a community).}}

\description{
This function provides the values of the second derivative of the beta profile. The second derivative represents the deceleration of the profile for any given value of "beta". The second derivative represents the deceleration of the profile for any given value of "beta". For "beta" close to ???1, high values of the second derivative denote high dominance because the profile tends to be extremely curved when some prevalent species are present in a community. For the interval [0,1] the second derivative tends always towards zero.
}
\author{
Fabrizio Maturo, Francesca Fortuna, Tonio Di Battista
}

\source{
ENVIRONMENTAL MONITORING THROUGH FUNCTIONAL BIODIVERSITY TOOLS - T. Di Battista, F. Fortuna, F. Maturo - Ecological Indicators (2016) 60, 237-247 - online (2015). DOI: 10.1016/j.ecolind.2015.05.056.

PARAMETRIC FUNCTIONAL ANALYSIS OF VARIANCE FOR FISH BIODIVERSITY ASSESSMENT - T. Di Battista, F. Fortuna, F. Maturo - Journal of Environmental Informatics. Special issue of Marine and Freshwater Environmental Engineering and Management, Journal of Environmental Informatics. In press.

}
\examples{
a=c(0.35,0.35,0.27,0.01,0.02)
b=c(0.54,0.20,0.17,0.06,0.03)
c=c(0.35,0.35,0.30,0,0)
d=c(0.51,0.31,0.07,0.10,0.01)
e=c(0.40,0.20,0.10,0.30,0)
x=matrix(rbind(a,b,c,d,e),5,5)
x
betasecond(x)

}
