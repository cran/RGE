\name{plot.RGE}
\alias{plot.RGE}
\title{Plots for object obtained by RGE}
\description{
This functions makes the plots of bayes estimate (mean) and the highest posterior density intervals at 95\%, of predicted value of genotype, his stability variance, and his bayesian yield stability index. 
}
\usage{
\method{plot}{RGE}(
    x,
    labelg = "Predicted value",
    labelsv = "Stability variance",
    labelby = "Bayesian yield stability index",
    margin = c(1, 0.8, 0, 0.8), ...)
}
\arguments{
  \item{x}{an object obtained by the function RGE \code{"RGE"}}
  \item{labelg}{Label to use in the plot of predicted value of genotype}
  \item{labelsv}{Label to use in the plot of stability variance}
  \item{labelby}{Label to use in the plot of bayesian yield stability index}
  \item{margin}{A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.}
  \item{\dots}{Further arguments to be passed}
}
\value{Plot of the predicted values, stability variance and bayesian yield stability index}
\author{Jose Miguel Cotes Torres \email{jmcotes@unal.edu.co}}

\seealso{\code{plot},\code{plot.mcmc},\code{par}}
\examples{
  data(m1)
  plot(m1)
}