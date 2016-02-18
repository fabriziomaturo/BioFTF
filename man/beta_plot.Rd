\name{beta_plot}
\alias{beta_plot}
\title{Plot the beta profile for each community.}
\usage{
beta_plot(x)
}
\arguments{\item{x}{A data matrix with the abundance of the species (the rows are the communities and the columns indicate the species living in a community).}}

\description{
This function provides a plot of the diversity profiles of the considered communities.
The domain "beta" represents the sensibility to changes in the presence of rare species.
It is a decreasing and concave upward curve showing different values of biodiversity for each value of beta. The most frequently used indices of biodiversity are special cases of this function. The diversity profile displays a complete picture of diversity. Its plot allows the analyst to rank different communities. Indeed, the higher curve highlights biological populations with higher diversity and vice versa. The plot tends towards a straight line, decreasing from ???1 to 1, in the case of maximum equitability with few species, and it becomes more curved with the prevalence of a few species over the others. Difficulties arise if we compare communities with intersecting profiles. In this case, it is impossible to determine which curve is the highest only using the diversity profile.

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
beta_plot(x)
}
