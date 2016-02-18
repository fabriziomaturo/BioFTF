\name{alltools}
\alias{alltools}
\title{Apply the whole analysis to the dataset.}
\usage{
alltools(x)
}
\arguments{\item{x}{A data matrix with the abundance of the species (the rows are the communities and the columns indicate the species living in a community).}}

\description{
This function provides the derivatives, the curvature, the radius of curvature and the arc length of diversity profiles. Morever, it computes the surface area under the beta profile and ranks the communities.
}

\author{
Fabrizio Maturo, Francesca Fortuna, Tonio Di Battista
}

\source{
ENVIRONMENTAL MONITORING THROUGH FUNCTIONAL BIODIVERSITY TOOLS - T. Di Battista, F. Fortuna, F. Maturo - Ecological Indicators (2016) 60, 237-247 - online (2015). DOI: 10.1016/j.ecolind.2015.05.056.

PARAMETRIC FUNCTIONAL ANALYSIS OF VARIANCE FOR FISH BIODIVERSITY ASSESSMENT - T. Di Battista, F. Fortuna, F. Maturo - Journal of Environmental Informatics. Special issue of Marine and Freshwater Environmental Engineering and Management, Journal of Environmental Informatics. In press.

}

\examples{
#example 1
x=matrix(c(3,5,5,2,1,2,6,8,5,1),2,5)
x
alltools(x)

#Some examples exspressed by relative abundance:

#example 2
x=matrix(c(0.3,0.5,0.1,0.05,0.05,0.25,0.25,0.25,0.25,0,0.35,0.3,0.35,0,0),3,5)
x
alltools(x)

#example 3
x=matrix(runif(1000),10,100)
#let's consider some missing species
x[1,1:20]=0
x[2,1:40]=0
x[3,1:60]=0
x[4,1:19]=0
x
alltools(x)

#example 4
x=matrix(runif(100),20,5)
x[1,1:2]=0
x[2,1:3]=0
x
alltools(x)

#example 5
a=c(0.35,0.35,0.27,0.01,0.02)
b=c(0.54,0.20,0.17,0.06,0.03)
c=c(0.35,0.35,0.30,0,0)
d=c(0.51,0.31,0.07,0.10,0.01)
e=c(0.40,0.20,0.10,0.30,0)
x=matrix(rbind(a,b,c,d,e),5,5)
x
alltools(x)
}

