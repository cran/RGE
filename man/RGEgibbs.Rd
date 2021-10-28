\name{RGEgibbs}
\alias{RGEgibbs}
\title{Response from Genotype to Environment}
\description{This function performance the GIBBS sampler for analyze reginals trials.
}
\usage{
RGEgibbs(data, gen_c, env_c, blk_c, y_c, prior.g = NULL,
        prior.vg =NULL,prior.b = NULL, prior.dfb = NULL,
        prior.sv = NULL, prior.dfsv = NULL, prior.se = NULL, 
        prior.dfse = NULL,
        burnin = 10, thin = 5, niter = 50, saveAt = 10)
}
\arguments{
  \item{data}{\code{data.frame}}
  \item{gen_c}{Number of the column from de data.frame with the genotypes information.}
  \item{env_c}{Number of the column from de data.frame with the environment information.}
  \item{blk_c}{Number of the column from de data.frame with the block information.}
  \item{y_c}{Number of the column from de data.frame with the phenotype information.}
  \item{prior.g}{Vector with prior information of the means of genotypes}
  \item{prior.vg}{Vector with prior information of the variances of the means of genotypes}
  \item{prior.b}{Vector with prior information of the variances of block within environment}
  \item{prior.dfb}{Vector with prior information of hyperparameter  degree of credibility of the variances of block within environment.}
  \item{prior.sv}{Vector with prior information of the stability's variances}
  \item{prior.dfsv}{Vector with prior information of hyperparameter degree of credibility of the stability's variances.}
  \item{prior.se}{Vector with prior information of the error's variances}
  \item{prior.dfse}{Vector with prior information of hyperparameter degree of credibility of the variances of error.}
  \item{burnin}{Number of iteration to be consider as burn-in period. This period is not saved in the final result.}
  \item{thin}{The thinning interval between consecutive observations. This interval is not saved in the final result.}
  \item{niter}{Numbers of iterations to be saved}
  \item{saveAt}{Save object outtS4 with samples of the posterior distribution on the work directory each "saveAt" iteration}
}
\value{Matrix with samples of the posterior distribution}
\references{
Cotes Torres, J. M., Gonzalez Jaimes, E. P., & Cotes Torres, A. (2016). Seleccion de Genotipos con Alta Respuesta y Estabilidad Fenotipica en Pruebas Regionales: Recuperando el Concepto Biologico. Revista Facultad De Ciencias Basicas, 8(2), 226-243.}
\author{Jose Miguel Cotes Torres \email{jmcotes@unal.edu.co}}
\examples{
##data(datos)
##m<-RGEgibbs(data=datos,gen_c=5,env_c=1,blk_c=4,y_c=9,
##thin=5,burnin=100,niter=10000,saveAt=1000)
}