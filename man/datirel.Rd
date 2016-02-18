\name{datirel}
\alias{datirel}
\title{Create a matrix with the relative frequencies for each community.}
\usage{datirel(x)}

\arguments{\item{x}{A data matrix with the absolute abundance of the species.}}

\description{
Create a matrix with the relative frequencies for each community.
}

\examples{
x=matrix(c(3,5,5,2,1,2,6,8,5,1),2,5)
x
datirel(x)

}
