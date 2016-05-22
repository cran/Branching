\name{rBGWM}
\alias{rBGWM}
\alias{rBGWM.gener}
\alias{rBGWM.indep}
\alias{rBGWM.multinom}
\title{Simulating a multi-type Bienayme - Galton - Watson process}
\description{Generate the trajectories of a multi-type Bienayme - 
Galton - Watson process from its offspring distributions, using three
different algorithms based on three different classes or families of these 
processes.}
\usage{
rBGWM(dists, type=c("general","multinomial","independents"), d,
      n, z0=rep(1,d), c.s=TRUE, tt.s=TRUE, rf.s=TRUE, file=NULL)
}
\arguments{
\item{dists }{offspring distributions. Its structure depends on the class
  of the Bienayme - Galton - Watson process (See details and examples).}
\item{type }{Class or family of the Bienayme - Galton - Watson process
  (See details).}
\item{d }{positive integer, number of types.}
\item{n }{positive integer, maximum lenght of the wanted trajectory.}
\item{z0 }{nonnegative integer vector of size d; initial population by type.}
\item{c.s }{logical value, if TRUE, the output object will include the 
generated trajectory of the process with the number of individuals for 
every combination parent type - descendent type.}
\item{tt.s }{logical value, if TRUE, the output object will include the 
generated trajectory of the process with the number of descendents by type.}
\item{rf.s }{logical value, if TRUE, the output object will include the 
generated trajectory of the process with the relative frequencies by type.}
\item{file }{the name of the output file where the generated trajectory of the 
process with the number of individuals for every combination parent type - 
descendent type could be stored.}
}
\details{

  This function performs a simulation of a multi-type Bienayme - Galton
  - Watson process (BGWM) from its offspring distributions.

  From particular offspring distributions and taking into account a
  differentiated algorithmic approach, we propose the following classes
  or types for these processes: 
  
  \code{general} This option is for BGWM processes without conditions over
  the offspring distributions, in this case, it is required as
  input data for each distribution, all d-dimensional vectors with their
  respective, greater than zero, probability.

  \code{multinomial} This option is for BGMW processes where each offspring
  distribution is a multinomial distribution with a random number of
  trials, in this case, it is required as input data, \eqn{d} univariate
  distributions related to the random number of trials for each
  multinomial distribution and a \eqn{d \times d} matrix where each row
  contains probabilities of the \eqn{d} possible outcomes for each multinomial
  distribution.
  
  \code{independents} This option is for BGMW processes where each offspring
  distribution is a joint distribution of \eqn{d} combined independent
  discrete random variables, one for each type of individuals, in this
  case, it is required as input data \eqn{d^2} univariate distributions.

  The structure need it for each classification is illustrated in the examples.

  These are the univariate distributions available:

  \emph{unif} Discrete uniform distribution, parameters \eqn{min} and
  \eqn{max}. All the non-negative integers between \eqn{min} y \eqn{max} have the same
  probability. 

  \emph{binom} Binomial distribution, parameters \eqn{n} and \eqn{p}.
  \deqn{p(x) = {n \choose x} {p}^{x} {(1-p)}^{n-x}}{
        p(x) = choose(n,x) p^x (1-p)^(n-x)}
  for x = 0, \eqn{\dots}{...}, n.

  \emph{hyper} Hypergeometric distribution, parameters \eqn{m} (the
  number of white balls in the urn), \eqn{n} (the number of white balls
  in the urn), \eqn{k} (the number of balls drawn from the urn). 
  \deqn{
    p(x) = \left. {m \choose x}{n \choose k-x} \right/ {m+n \choose k}%
  }{p(x) =      choose(m, x) choose(n, k-x) / choose(m+n, k)}
  for x = 0, ..., k.
  
  \emph{geom} Geometric distribution, parameter \eqn{p}.
  \deqn{
    p(x) = p {(1-p)}^{x}
  }{p(x) = p (1-p)^x}
  for x = 0, 1, 2, \eqn{\dots}{...}

  \emph{nbinom} Negative binomial distribution, parameters \eqn{n} and
  \eqn{p}.
  \deqn{
    p(x) = \frac{\Gamma(x+n)}{\Gamma(n) x!} p^n (1-p)^x
  }{p(x) = Gamma(x+n)/(Gamma(n) x!) p^n (1-p)^x}
  for x = 0, 1, 2, \eqn{\dots}{...}

  \emph{pois} Poisson distribution, parameter \eqn{\lambda}{lambda}.
  \deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{x!}}{p(x) = lambda^x exp(-lambda)/x!}
  for x = 0, 1, 2, \eqn{\dots}{...}

  \emph{norm} Normal distribution rounded to integer values and negative
  values become 0, parameters \eqn{\mu}{mu} and
  \eqn{\sigma}{sigma}.
  \deqn{
    p(x) = \int_{x-0.5}^{x+0.5} \frac{1}{\sqrt{2\pi}\sigma} e^{-(t-\mu)^2/2\sigma^2}dt%
  }{p(x) = \int_{x-0.5}^{x+0.5} 1/(sqrt(2 pi) sigma) e^-((t - mu)^2/(2sigma^2)) dt}
  for x = 1, 2, \eqn{\dots}{...}
  \deqn{p(x) = \int_{-\infty}^{0.5} \frac{1}{\sqrt{2\pi}\sigma} e^{-(t-\mu)^2/2\sigma^2}dt%
  }{p(x) = \int_{-\infty}^{0.5} 1/(sqrt(2 pi) sigma) e^-((t - mu)^2/(2sigma^2)) dt}
  for x = 0

  \emph{lnorm} Lognormal distribution rounded to integer values,
  parameters \code{logmean} \eqn{=\mu}{= mu} y \code{logsd} \eqn{=\sigma}{= sigma}.
  \deqn{
    p(x) = \int_{x-0.5}^{x+0.5} \frac{1}{\sqrt{2\pi}\sigma t} e^{-(\log(t) - \mu)^2/2 \sigma^2 }dt%
  }{p(x) = \int_{x-0.5}^{x+0.5} 1/(sqrt(2 pi) sigma t) e^-((log t - mu)^2 / (2sigma^2))dt}
  for x = 1, 2, \eqn{\dots}{...}
  \deqn{
    p(x) = \int_{0}^{0.5} \frac{1}{\sqrt{2\pi}\sigma t} e^{-(\log(t) - \mu)^2/2 \sigma^2 }dt%
  }{p(x) = \int_{0}^{0.5} 1/(sqrt(2 pi) sigma t) e^-((log t - mu)^2 / (2sigma^2)) dt}
  for x = 0

  \emph{gamma} Gamma distribution rounded to integer values,
  parameters \code{shape} \eqn{=\alpha}{= a} y \code{scale} \eqn{=\sigma}{= s}.
  \deqn{
    p(x)= \int_{x-0.5}^{x+0.5} \frac{1}{{\sigma}^{\alpha}\Gamma(\alpha)}{t}^{\alpha-1} e^{-t/\sigma} dt%
  }{p(x)= \int_{x-0.5}^{x+0.5} 1/(s^a Gamma(a)) t^(a-1) e^-(t/s) dt}
  para x = 1, 2, \eqn{\dots}{...}
  \deqn{
    p(x)= \int_{0}^{0.5} \frac{1}{{\sigma}^{\alpha}\Gamma(\alpha)}{t}^{\alpha-1} e^{-t/\sigma} dt%
  }{p(x)= \int_{0}^{0.5} 1/(s^a Gamma(a)) t^(a-1) e^-(t/s) dt}
  for x = 0
}
\value{
An object of class \code{list} with these components:
\item{i.d }{input. number of types.}
\item{i.dists }{input. offspring distributions.}
\item{i.n }{input. maximum lenght of the generated trajectory.}
\item{i.z0 }{input. initial population by type.}
\item{o.c.s }{output. A matrix indicating the generated trajectory of the 
process with the number of individuals for every combination parent type - 
descendent type.}
\item{o.tt.s}{output. A matrix indicating the generated trajectory of the 
process with the number of descendents by type.}
\item{o.rf.s}{output. A matrix indicating the generated trajectory of the 
process with the relative frequencies by type.}
}
\references{
Torres-Jimenez, C. J. (2010), Relative frequencies and parameter 
estimation in multi-type Bienaym? - Galton - Watson processes, 
Master's Thesis, Master of Science in Statistics.
Universidad Nacional de Colombia. Bogota, Colombia.

Stefanescu, C. (1998), 'Simulation of a multitype Galton-Watson
chain', Simulation Practice and Theory 6(7), 657-663.

Athreya, K. & Ney, P. (1972), Branching Processes, Springer-Verlag.
}
\author{
Camilo Jose Torres-Jimenez \email{cjtorresj@unal.edu.co}
}
\seealso{
\code{\link{BGWM.mean}}, \code{\link{BGWM.covar}}, \code{\link{BGWM.mean.estim}}, \code{\link{BGWM.covar.estim}}
}
\examples{
\dontrun{
## Simulation based on a model analyzed in Stefanescu(1998)

# Variables and parameters
d <- 2
n <- 30
N <- c(90, 10)
a <- c(0.2, 0.3)

# with independent distributions
Dists.i <- data.frame( name=rep( "pois", d*d ),
                       param1=rep( a, rep(d,d) ),
                       stringsAsFactors=FALSE )

rA <- rBGWM(Dists.i, "independents", d, n, N)

# with multinomial distributions
dist <- data.frame( name=rep( "pois", d ),
                    param1=a*d,
                    stringsAsFactors=FALSE )
matrix.b <- matrix( rep(0.5, 4), nrow=2 )
Dists.m <- list( dists.eta=dist, matrix.B=matrix.b )

rB <- rBGWM(Dists.m, "multinomial", d, n, N)

# with general distributions (approximation)
max <- 30
A <- t(expand.grid(c(0:max),c(0:max)))
aux1 <- factorial(A)
aux1 <- apply(aux1,2,prod)
aux2 <- apply(A,2,sum)
distp <- function(x,y,z){ exp(-d*x)*(x^y)/z }
p <- sapply( a, distp, aux2, aux1 )
prob <- list( dist1=p[,1], dist2=p[,2] )
size <- list( dist1=ncol(A), dist2=ncol(A) )
vect <- list( dist1=t(A), dist2=t(A) )
Dists.g <- list( sizes=size, probs=prob, vects=vect )

rC <- rBGWM(Dists.g, "general", d, n, N)

# Comparison chart
dev.new()
plot.ts(rA$o.tt.s,main="with independents")
dev.new()
plot.ts(rB$o.tt.s,main="with multinomial")
dev.new()
plot.ts(rC$o.tt.s,main="with general (aprox.)")
}
}
\keyword{Bienayme - Galton - Watson}
\keyword{multytipe}
\keyword{branching}
\keyword{generate}
\keyword{random}
