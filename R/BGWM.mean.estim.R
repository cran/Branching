# These functions calculate a mean estimation from observed sample of 
# Bienayme - Galton - Watson multitype processes. 
# Copyright (C) 2010  Camilo Jose Torres-Jimenez <cjtorresj@unal.edu.co>

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

BGWM.mean.estim <- function(sample, method=c("EE","MLE"), d, n, z0)
{
  
  method <- match.arg(method)
  method <- switch(method,
                   "EE" = 1,
                   "MLE" = 2)

  m <- switch(method,
		 {#1
		   m <- BGWM.mean.EE(sample, d, n, z0)
		   m
		 },
		 {#2
		   m <- BGWM.mean.MLE(sample, d, n, z0)
		   m
		 })
  
  colnames(m) <- paste( "type", 1:d, sep="" )
  rownames(m) <- paste( "type", 1:d, sep="" )

  list( method=switch( method, "Empirical Estimation of the means", "Maximum Likelihood Estimation of the means" ), m=m )
  
}  


BGWM.mean.EE <- function(y, d, n, z0)
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

  z <- y[seq( (n-2)*d+1 , (n-1)*d, 1 ),]
  if(n != 1)
    z <- apply( z, 2, sum )
  else
    z <- z0
  z <- diag( (1 / z) )
  y <- y[seq( (n-1)*d+1 , n*d, 1 ),]
  out <- z %*% y

  out

}



BGWM.mean.MLE <- function(y, d, n, z0)
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

  if(n != 1)
  {
    z <- rbind( z0, y[seq( 1, (n-1)*d, 1 ),] )
    z <- apply( z, 2, sum )
  }
  else
    z <- z0
  z <- diag( ( 1 / z ) )
  y <- as.data.frame(y[seq( 1, (n*d), 1),])
  y[,"type"] <- as.factor( rep( 1:d, n ) )
  y <- aggregate( y[, 1:d], list( y[,"type"] ), sum )[,-1]
  out <- as.matrix(z) %*% as.matrix(y)
  
  out

}