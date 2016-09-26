\name{area}
\alias{area}
\title{Compute the area under the beta profile for each community.}
\usage{
area(x,n)
}

\arguments{

\item{x}{A data matrix with the abundance of the species (the rows are the communities and the columns indicate the species living in a community). The user can freely choose to use the absolute or relative frequencies.}

\item{n}{The number of points of the domain [-1,1] that the user wants to consider in computing the functional tools. It is suggested to choose a number in the interval [3,10000] because a value of less than 3 has little sense; indeed, the major indices of biodiversity are the richness index, the index of Shannon and the index of Simpson (that we get for beta equal to -1, 0, +1, respectively). On the other hand, a value greater than 10000 is not very interesting because with a value of approximately 100 we already have a good approximation. The examples provided in this package always divide the domain into 20 intervals. The default value is n=20.}



}


\description{
This function provides the area under the beta profile for each community. In an ecological framework, the area under the profiles can be used to assess biodiversity because it can rank communities, even if they have different number of species. Indeed, in a case of maximum dominance, the area is extremely low; while in a case of evenness, it increase. The main advantage of using area is that the ordering among communities can be investigated without the analysis of a graph. It provides a scalar measure of biodiversity considering the whole domain. This is not possible with the classical indices.
}

\examples{
x=matrix(c(0.3,0.5,0.1,0.05,0.05,0.25,0.25,0.25,0.25,0,0.35,0.3,0.35,0,0),3,5)
area(x,20)

}
