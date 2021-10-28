\name{summary.RGE}
\alias{summary.RGE}
\title{Summary of the posterior distribution}
\description{
Obtaining the bayes estimative and the highest posterior density intervals at 95\% from predicted value of the genotypes, his stability variances and his bayesian yiled stability indexes.
}
\usage{
\method{summary}{RGE}(object, ...)
}
\arguments{
  \item{object}{An object obtained by the function RGE \code{"RGE"}}
  \item{\dots}{Further arguments to be passed}
}
\value{Return a list with:
  \item{mu}{Summary with predicted values of genotypes.}
  \item{sv}{Summary of the stability variances}
  \item{sv}{Summary of the bayesian yield stability indexes}
}
\author{Jose Miguel Cotes Torres \email{jmcotes@unal.edu.co}}
\note{coda package is needed.}
\seealso{\code{coda}}
\examples{
  data(m1)
  summary(m1)
}