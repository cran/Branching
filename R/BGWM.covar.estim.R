# These functions calculate a covariance estimation from observed sample of 
# Bienayme - Galton - Watson multitype processes. 
# Copyright (C) 2010  Camilo José Torres Jiménez <cjtorresj@unal.edu.co>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


BGWM.covar.estim <- function(sample, method=c("EE-m","MLE-m"), d, n, z0)
{
  
  method <- match.arg(method)
  method <- switch(method,
                   "EE-m" = 1,
                   "MLE-m" = 2)

  V <- switch(method,
		 {#1
		   V <- BGWM.covar.EE(sample, d, n, z0)
		   V
		 },
		 {#2
		   V <- BGWM.covar.MLE(sample, d, n, z0)
		   V
		 })
  
  dimnames(V) <- list( paste( "dist", rep(1:d,rep(d,d)), ".type", rep(1:d,d), sep="" ), paste( "type", 1:d, sep="" ) )

  list(method=switch( method, "with Empirical Estimation of the means", "with Maximum Likelihood Estimation of the means" ), V=V )
  
}  


BGWM.covar.EE <- function(y, d, n, z0)
{

  y <- as.matrix(y) 
  if(length(d) != 1)
    stop("'d' must be a number")
  if(length(n) != 1)
    stop("'n' must be a number")
  if(length(z0) != d)
    stop("'z0' must be a d-dimensional vector")
  if(TRUE %in% (z0 < 0))
    stop("'z0' must have positive elements")
  if(is.matrix(y) == FALSE)
    stop("'y' must be a matrix")
  if(ncol(y) != d || nrow(y) < (n*d))
    stop("'y' must have d columns and at least (n*d) rows")
  if(n == 1)
    stop("'n' must be greater than 1")

  y <- y[1:(n*d),]
  out <- matrix( rep( 0, (d*d*d) ), ncol=d )
  Mn <- BGWM.mean.EE(y, d, n, z0)

  for( i in 1:(n-1) )
  {
    if(i != 1)
      Zi_1 <- apply( y[seq( (i-2)*d+1, (i-1)*d, 1 ),], 2, sum )
    else
      Zi_1 <- z0
    Zi_1 <- rep( Zi_1, rep( d, d ) )
    aux <- BGWM.mean.EE(y, d, i, z0) - Mn
    out <- out + matrix( c( apply( aux, 1, tcrossprod ) ), ncol=d, byrow=TRUE ) * Zi_1
  }

  out <- out / n

  out

}


BGWM.covar.MLE <- function(y, d, n, z0)
{

  y <- as.matrix(y) 
  if(length(d) != 1)
    stop("'d' must be a number")
  if(length(n) != 1)
    stop("'n' must be a number")
  if(length(z0) != d)
    stop("'z0' must be a d-dimensional vector")
  if(TRUE %in% (z0 < 0))
    stop("'z0' must have positive elements")
  if(is.matrix(y) == FALSE)
    stop("'y' must be a matrix")
  if(ncol(y) != d || nrow(y) < (n*d))
    stop("'y' must have d columns and at least (n*d) rows")
  if(n == 1)
    stop("'n' must be greater than 1")

  y <- y[1:(n*d),]

  out <- matrix( rep( 0, (d*d*d) ), ncol=d )
  Mn <- BGWM.mean.MLE(y, d, n, z0)

  for( i in 1:(n-1) )
  {
    if(i != 1)
      aux2 <- aux2 + apply( y[seq( (i-2)*d+1, (i-1)*d, 1 ),], 2, sum )
    else
      aux2 <- z0
    aux1 <- BGWM.mean.MLE(y, d, i, z0) - Mn
    out <- out + matrix( c( apply( aux1, 1, tcrossprod ) ), ncol=d, byrow=TRUE ) * rep( aux2, rep( d, d ) )
  }

  out <- out / n

  out

}