\name{BGWM.mean.estim}
\alias{BGWM.mean.estim}
\alias{BGWM.mean.EE}
\alias{BGWM.mean.MLE}
\title{Estimation of the mean matrix of a multi-type Bienayme
  - Galton - Watson process} 
\description{Calculates a estimation of the mean matrix of a multi-type Bienayme
  - Galton - Watson process from experimental observed data that can be
  modeled by this kind of process.}
\usage{
BGWM.mean.estim(sample, method=c("EE","MLE"), d, n, z0)
}
\arguments{
\item{sample }{nonnegative integer matrix with \eqn{d} columns and
  \eqn{dn} rows, trajectory
  of the process with the number of individuals for every 
  combination parent type - descendent type (observed data).}
\item{method }{methods of estimation (EE Empirical estimacion,
  MLE Maximum likelihood estimation).}
\item{d }{positive integer, number of types.}
\item{n }{positive integer, nth generation.}
\item{z0 }{nonnegative integer vector of size d, initial population by
  type.}
}
\details{
  This function estimates the mean matrix of a BGWM process using two
  possible estimators, empirical estimator and maximum likelihood
  estimator, they both require the so-called full sample associated with
  the process, ie, it is required to have the trajectory of the process
  with the number of individuals for every combination parent type -
  descendent type. For more details see Torres-Jimenez (2010) or
  Maaouia & Touati (2005).
}
\value{
  A \code{list} object with:
  \item{method }{method of estimation selected.}
  \item{m }{A \code{matrix} object, estimation of the \eqn{d \times d}
    mean matrix of the process.}
}
\references{
Torres-Jimenez, C. J. (2010), Relative frequencies and parameter 
estimation in multi-type Bienayme - Galton - Watson processes, 
Master's Thesis, Master of Science in Statistics.
Universidad Nacional de Colombia. Bogota, Colombia.

Maaouia, F. & Touati, A. (2005), 'Identification of Multitype Branching
Processes', The Annals of Statistics 33(6), 2655-2694.
}
\author{
Camilo Jose Torres-Jimenez \email{cjtorresj@unal.edu.co}
}
\seealso{
\code{\link{BGWM.mean}}, \code{\link{BGWM.covar}}, \code{\link{rBGWM}}, \code{\link{BGWM.covar.estim}}
}
\examples{
\dontrun{
## Estimation of mean matrix from simulated data

# Variables and parameters
d <- 3
n <- 30
N <- c(10,10,10)
LeslieMatrix <- matrix( c(0.08, 1.06, 0.07, 
                          0.99, 0, 0, 
                          0, 0.98, 0), 3, 3 )

# offspring distributions from the Leslie matrix
# (with independent distributions)  
Dists.pois <- data.frame( name=rep( "pois", d ),
                          param1=LeslieMatrix[,1],
                          param2=NA,
                          stringsAsFactors=FALSE )
Dists.binom <- data.frame( name=rep( "binom", 2*d ),
                           param1=rep( 1, 2*d ),
                           param2=c(t(LeslieMatrix[,-1])),
                           stringsAsFactors=FALSE ) 
Dists.i <- rbind(Dists.pois,Dists.binom)
Dists.i <- Dists.i[c(1,4,5,2,6,7,3,8,9),]
Dists.i

# mean matrix of the process from its offspring distributions
m <- BGWM.mean(Dists.i,"independents",d)

# generated trajectories of the process from its offspring distributions
simulated.data <- rBGWM(Dists.i, "independents", d, n, N, 
                        TRUE, FALSE, FALSE)$o.c.s

# mean matrix empiric estimate from generated trajectories of the process
m.EE <- BGWM.mean.estim( simulated.data, "EE", d, n, N )$m

# mean matrix maximum likelihood estimate from generated trajectories
# of the process 
m.MLE <- BGWM.mean.estim( simulated.data, "MLE", d, n, N )$m

# Comparison of exact and estimated mean matrices
m
m - m.EE
m - m.MLE
}
}
\keyword{Bienayme - Galton - Watson}
\keyword{multi-type}
\keyword{branching}
\keyword{estimation}
\keyword{mean}