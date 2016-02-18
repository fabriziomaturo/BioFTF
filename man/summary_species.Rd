\name{summary_species}
\alias{summary_species}
\title{Create a species summary for each community. Based on the dataset with absolute abundance}
\usage{summary_species(x)}

\arguments{\item{x}{A data matrix with the absolute abundance of the species.}}

\description{
Create the species summary for each community (absolute abundance).
}

\examples{

x=matrix(c(3,5,5,2,1,2,6,8,5,1),2,5)
x
summary_species(x)


}
