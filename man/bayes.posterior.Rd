\name{bayes.posterior}
\alias{bayes.posterior}
\title{Summary of the posterior distribution}
\description{
Obtaining the bayes estimative and the highest posterior density intervals at 95\% from a object obtained by RGE function.
}
\usage{
bayes.posterior(x, ...)
}
\arguments{
  \item{x}{An object obtained by the function RGE \code{"RGE"}}
  \item{\dots}{Further arguments to be passed}
}
\value{Dataframe with the summary of the posterior distribution}
\author{
Jose Miguel Cotes Torres \email{jmcotes@unal.edu.co}
}
\note{coda package is needed.}
\seealso{\code{coda}}
\examples{
data(m1)
bayes.posterior(m1)
}