\name{BGWM.covar.estim}
\alias{BGWM.covar.estim}
\alias{BGWM.covar.EE}
\alias{BGWM.covar.MLE}
\title{Estimation of the covariance matrices of a multi-type Bienayme
  - Galton - Watson process}
\description{Calculates a estimation of the covariance matrices of a
  multi-type Bienayme - Galton - Watson process from experimental
  observed data that can be modeled by this kind of process.} 
\usage{
BGWM.covar.estim(sample, method=c("EE-m","MLE-m"), d, n, z0)
}
\arguments{
\item{sample }{nonnegative integer matrix with \eqn{d} columns and
  \eqn{dn} rows, trajectory
  of the process with the number of individuals for every 
  combination parent type - descendent type (observed data).}
\item{method }{methods of estimation (EE-m with empirical estimation
  of the mean matrix, MLE-m with maximum likelihood estimation of the
  mean matrix).}
\item{d }{positive integer, number of types.}
\item{n }{positive integer, nth generation.}
\item{z0 }{nonnegative integer vector of size d, initial population by
  type.}
}
\details{
  This function estimates the covariance matrices of a BGWM process using two
  possible estimators from asymptotic results related with empirical
  estimator and maximum likelihood estimator of the mean matrix, they
  both require the so-called full sample associated with the process,
  ie, it is required to have the trajectory of the process with the
  number of individuals for every combination parent type - descendent
  type. For more details see Torres-Jimenez (2010) or Maaouia & Touati
  (2005). 
}
\value{
  A \code{list} object with:
  \item{method }{method of estimation selected.}
  \item{V }{A \code{matrix} object, estimation of the \eqn{d} covariance
    matrices of the process, combined by rows.}
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
\code{\link{BGWM.mean}}, \code{\link{BGWM.covar}}, \code{\link{BGWM.mean.estim}}, \code{\link{rBGWM}}
}
\examples{
\dontrun{
## Estimation of covariace matrices from simulated data

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

# covariance matrices of the process from its offspring distributions
V <- BGWM.covar(Dists.i,"independents",d)

# generated trajectories of the process from its offspring distributions
simulated.data <- rBGWM(Dists.i, "independents", d, n, N, 
                        TRUE, FALSE, FALSE)$o.c.s

# estimation of covariance matrices using mean matrix empiric estimate
# from generated trajectories of the process 
V.EE.m <- BGWM.covar.estim( simulated.data, "EE-m", d, n, N )$V

# estimation of covariance matrices using mean matrix maximum likelihood
# estimate from generated trajectories of the process 
V.MLE.m <- BGWM.covar.estim( simulated.data, "MLE-m", d, n, N )$V

# Comparison of exact and estimated covariance matrices
V
V - V.EE.m
V - V.MLE.m
}
}
\keyword{Bienayme - Galton - Watson}
\keyword{multi-type}
\keyword{branching}
\keyword{estimation}
\keyword{variance}
\keyword{covariance}