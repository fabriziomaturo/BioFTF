\name{summary_species_relative}
\alias{summary_species_relative}
\title{Create a species summary for each community. Based on the dataset with relative abundance. }
\usage{summary_species_relative(x)}

\arguments{\item{x}{A data matrix with the absolute abundance of the species.}}

\description{
Create a species summary for each community (relative abundance).
}

\examples{

x=matrix(c(3,5,5,2,1,2,6,8,5,1),2,5)
x
summary_species_relative(x)

}
