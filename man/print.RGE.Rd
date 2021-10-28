\name{print.RGE}
\alias{print.RGE}
\title{Summary of the posterior distribution}
\description{
Obtaining the bayes estimative and the highest posterior density intervals at 95\% from predicted value of the genotypes, his stability variances and his bayesian yiled stability indexes.
}
\usage{
\method{print}{RGE}(x, ...)
}
\arguments{
  \item{x}{An object obtained by the function RGE \code{"RGE"}}
  \item{\dots}{Further arguments to be passed}
}
\value{Do not return any value. It is a print version of summary.RGE}
\author{Jose Miguel Cotes Torres \email{jmcotes@unal.edu.co}}
\note{coda package is needed.}
\seealso{\code{coda}}
\examples{
  data(m1)
  print(m1)
}