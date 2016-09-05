\name{beta}
\alias{beta}
\title{Compute the values of beta profile for each community in every part of the beta domain.}
\usage{
beta(x,n)
}


\arguments{

\item{x}{A data matrix with the abundance of the species (the rows are the communities and the columns indicate the species living in a community). The user can freely choose to use the absolute or relative frequencies.}

\item{n}{The number of points of the domain [-1,1] that the user wants to consider in computing the functional tools. It is suggested to choose a number in the interval [3,10000] because a value of less than 3 has little sense; indeed, the major indices of biodiversity are the richness index, the index of Shannon and the index of Simpson (that we get for beta equal to -1, 0, +1, respectively). On the other hand, a value greater than 10000 is not very interesting because with a value of approximately 100 we already have a good approximation. The examples provided in this package always divide the domain into 20 intervals. The default value is n=20, thus it is also possible to omit the value of n.}



}


\description{
The matrix obtained with "beta" contains the values of diversity measures computed for each community by column and for each point of the domain. This matrix provides the values of diversity profile; it displays a complete picture of diversity. Indeed, the most frequently used indices of biodiversity are special cases of this function.


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
beta(x,20)

}
