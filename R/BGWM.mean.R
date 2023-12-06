# These functions calculate the mean of a Bienayme - Galton - Watson 
# multitype processes. 
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


BGWM.mean <- function(dists, type=c("general","multinomial","independents"), d, n=1, z0=NULL, maxiter = 1e5)
{
  
  if(n < 1)
    stop("'n' must be a positive number")
  #parameters restrictions
  
  type <- match.arg(type)
  type <- switch(type,
                 "general" = 1,
                 "multinomial" = 2,
                 "independents" = 3)

  M <- switch(type,
		 {#1
		   M <- BGWM.gener.mean(dists, d)
		   M
		 },
		 {#2
		   M <- BGWM.multinom.mean(dists, d, maxiter)
		   M
		 },
		 {#3
		   M <- BGWM.indep.mean(dists, d, maxiter)
		   M
		 })
  
  dimnames(M) <- list( paste( "type", 1:d, sep="" ), paste( "type", 1:d, sep="" ) )

  mean <- M
  if( n != 1 )
  {
    if( n > 1 )
      # This could be changed to a better method  
      for( i in 2:n )
	mean <- mean %*% M
  }

  if( length(z0) != 0 )
    mean <- z0 %*% mean

  mean

}



BGWM.gener.mean <- function(gener.dists, d)
{

  #parameters restrictions
 
  s <- gener.dists[[1]]
  p <- gener.dists[[2]]
  v <- gener.dists[[3]]

  probs <- data.frame( p = unlist(p), k = as.factor( rep( 1:d, unlist(s) ) ) )
  vectors <- data.frame( v = matrix( unlist( lapply( v, t ) ), ncol=d, byrow=TRUE), k = as.factor( rep( 1:d, unlist(s) ) ) )
  probs[,1] <- probs[,1] / rep( aggregate(probs[,1], list(probs[,2]) , sum)[,2] , unlist(s) )
  vectors[,-(d+1)] <- vectors[,-(d+1)] * probs[,1]
  mean <- as.matrix(aggregate( vectors[,-(d+1)], list(vectors[,(d+1)]), sum )[,-1])

  mean

}


BGWM.multinom.mean <- function(multinom.dists, d, maxiter = 1e5)
{

  #parameters restrictions

  dists <- multinom.dists[[1]]
  pmultinom <- multinom.dists[[2]]
  pmultinom <- pmultinom / apply(pmultinom, 1, sum)
  mean <- rep( NA, d )

  # unif
  a <- dists[,1] == "unif"
  if(TRUE %in% a)
    mean[a] <- ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) ) / 2
  # binom
  a <- dists[,1] == "binom"
  if(TRUE %in% a)
    mean[a] <- as.numeric( dists[a,2] ) * as.numeric( dists[a,3] )
  # hyper
  a <- dists[,1] == "hyper"
  if(TRUE %in% a)
    mean[a] <- as.numeric( dists[a,2] ) * as.numeric( dists[a,4] ) / ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) )
  # geom
  a <- dists[,1] == "geom"
  if(TRUE %in% a)
    mean[a] <- ( 1 - as.numeric( dists[a,2] ) ) / as.numeric( dists[a,2] )
  # nbinom
  a <- dists[,1] == "nbinom"
  if(TRUE %in% a)
    mean[a] <- as.numeric( dists[a,2] ) * ( 1 - as.numeric( dists[a,3] ) ) / as.numeric( dists[a,3] )
  # pois
  a <- dists[,1] == "pois"
  if(TRUE %in% a)
    mean[a] <- as.numeric( dists[a,2] )
  # norm
  a <- dists[,1] == "norm"
  n <- length( mean[a] )
  if(n > 0)
    mean[a] <- round( .C("param_estim_roundcut0_norm",
			  as.integer( maxiter ),
			  as.integer( n ),
			  as.double( as.numeric( dists[a,2] ) ),
			  as.double( as.numeric( dists[a,3] ) ),
			  as.integer( 1 ),
			  mean.estim=double( n ),
			  var.estim=double( n ),
			  PACKAGE="Branching")$mean.estim , round( log10(maxiter)/2 ) )
  # lnorm
  a <- dists[,1] == "lnorm"
  n <- length( mean[a] )
  if(n > 0)
    mean[a] <- round( .C("param_estim_round_lnorm",
			  as.integer( maxiter ),
			  as.integer( n ),
			  as.double( as.numeric( dists[a,2] ) ),
			  as.double( as.numeric( dists[a,3] ) ),
			  as.integer( 1 ),
			  mean.estim=double( n ),
			  var.estim=double( n ),
			  PACKAGE="Branching")$mean.estim , round( log10(maxiter)/2 ) )
  # gamma
  a <- dists[,1] == "gamma"
  n <- length( mean[a] )
  if(n > 0)
    mean[a] <- round( .C("param_estim_round_gamma",
			  as.integer( maxiter ),
			  as.integer( n ),
			  as.double( as.numeric( dists[a,2] ) ),
			  as.double( as.numeric( dists[a,3] ) ),
			  as.integer( 1 ),
			  mean.estim=double( n ),
			  var.estim=double( n ),
			  PACKAGE="Branching")$mean.estim , round( log10(maxiter)/2 ) )

  mean <- diag( mean ) %*% pmultinom

  mean

}


BGWM.indep.mean <- function(indep.dists, d, maxiter = 1e5)
{
  
  #parameters restrictions

  dists <- indep.dists 
  mean <- rep( NA, (d*d) )

  # unif
  a <- dists[,1] == "unif"
  if(TRUE %in% a)
    mean[a] <- ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) ) / 2
  # binom
  a <- dists[,1] == "binom"
  if(TRUE %in% a)
    mean[a] <- as.numeric( dists[a,2] ) * as.numeric( dists[a,3] )
  # hyper
  a <- dists[,1] == "hyper"
  if(TRUE %in% a)
    mean[a] <- as.numeric( dists[a,2] ) * as.numeric( dists[a,4] ) / ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) )
  # geom
  a <- dists[,1] == "geom"
  if(TRUE %in% a)
    mean[a] <- ( 1 - as.numeric( dists[a,2] ) ) / as.numeric( dists[a,2] )
  # nbinom
  a <- dists[,1] == "nbinom"
  if(TRUE %in% a)
    mean[a] <- as.numeric( dists[a,2] ) * ( 1 - as.numeric( dists[a,3] ) ) / as.numeric( dists[a,3] )
  # pois
  a <- dists[,1] == "pois"
  if(TRUE %in% a)
    mean[a] <- as.numeric( dists[a,2] )
  # norm
  a <- dists[,1] == "norm"
  n <- length( mean[a] )
  if(n > 0)
    mean[a] <- round( .C("param_estim_roundcut0_norm",
			  as.integer( maxiter ),
			  as.integer( n ),
			  as.double( as.numeric( dists[a,2] ) ),
			  as.double( as.numeric( dists[a,3] ) ),
			  as.integer( 1 ),
			  mean.estim=double( n ),
			  var.estim=double( n ),
			  PACKAGE="Branching")$mean.estim , round( log10(maxiter)/2 ) )
  # lnorm
  a <- dists[,1] == "lnorm"
  n <- length( mean[a] )
  if(n > 0)
    mean[a] <- round( .C("param_estim_round_lnorm",
			  as.integer( maxiter ),
			  as.integer( n ),
			  as.double( as.numeric( dists[a,2] ) ),
			  as.double( as.numeric( dists[a,3] ) ),
			  as.integer( 1 ),
			  mean.estim=double( n ),
			  var.estim=double( n ),
			  PACKAGE="Branching")$mean.estim , round( log10(maxiter)/2 ) )
  # gamma
  a <- dists[,1] == "gamma"
  n <- length( mean[a] )
  if(n > 0)
    mean[a] <- round( .C("param_estim_round_gamma",
			  as.integer( maxiter ),
			  as.integer( n ),
			  as.double( as.numeric( dists[a,2] ) ),
			  as.double( as.numeric( dists[a,3] ) ),
			  as.integer( 1 ),
			  mean.estim=double( n ),
			  var.estim=double( n ),
			  PACKAGE="Branching")$mean.estim , round( log10(maxiter)/2 ) )

  mean <- matrix( mean, nrow=d, byrow=TRUE )

  mean

}