\name{betaprime}
\alias{betaprime}
\title{Compute the values of the first derivative of beta profile for each community in every part of the beta domain.}
\usage{
betaprime(x,n)
}

\arguments{

\item{x}{A data matrix with the abundance of the species (the rows are the communities and the columns indicate the species living in a community). The user can freely choose to use the absolute or relative frequencies.}

\item{n}{The number of points of the domain [-1,1] that the user wants to consider in computing the functional tools. It is suggested to choose a number in the interval [3,10000] because a value of less than 3 has little sense; indeed, the major indices of biodiversity are the richness index, the index of Shannon and the index of Simpson (that we get for beta equal to -1, 0, +1, respectively). On the other hand, a value greater than 10000 is not very interesting because with a value of approximately 100 we already have a good approximation. The examples provided in this package always divide the domain into 20 intervals. The default value is n=20, thus it is also possible to omit the value of n.}



}


\description{
This function provides the values of the first derivative of beta profile. This is an essential component of functional data analysis (Ramsay and Silverman, 2005), and in an ecological framework, is useful for highlighting the characteristics of the profiles and explaining some of the variation from curve to curve.
For "beta" close to -1, high absolute values for the first derivative correspond to a strong decrease of the profile; thus, a sudden decrease suggests the presence of few prevalent species in a community (dominance). On the contrary, low absolute values indicate that living organisms tend to be equally distributed among the species (evenness).
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
betaprime(x,20)

}
