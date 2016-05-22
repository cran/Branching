# These functions calculate the covariance from input variables of 
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


BGWM.covar <- function(dists, type=c("general","multinomial","independents"), d, n=1, z0=NULL, maxiter = 1e5)
{
  
  if(n < 1)
    stop("'n' must be a positive number")
  #parameters restrictions

  stype <- match.arg(type)
  stype <- switch(stype,
                  "general" = 1,
                  "multinomial" = 2,
                  "independents" = 3)
  
  V <- switch( stype,
		 {#1
		   V <- BGWM.gener.covar(dists, d)
		   V
		 },
		 {#2
		   V <- BGWM.multinom.covar(dists, d, maxiter)
		   V
		 },
		 {#3
		   V <- BGWM.indep.covar(dists, d, maxiter)
		   V
		 })
  
  if( n > 1 )
  {
    m <- BGWM.mean( dists, type, d, maxiter=maxiter )  
    V <- matrix( t(V), d*d, d )
    m.n_i <- diag( rep( 1, d ), d, d )
    if( length(z0) != 0 )
    {
      Cov <- rowSums( V %*% diag( c( BGWM.mean( dists, type, d, (n-1), z0, maxiter ) ), d, d ) )
      for(i in (n-1):1)
      {
	m.n_i <- m.n_i %*% m
	if( i != 1 )
	  AUX <- rowSums( V %*% diag( c( BGWM.mean( dists, type, d, (i-1), z0, maxiter ) ), d, d ) )
	else
	  AUX <- rowSums( V %*% diag( z0, d, d ) )
	AUX <- matrix( AUX, d, d )
	AUX <- t(m.n_i) %*% AUX %*% m.n_i
	Cov <- Cov + AUX
      }
      Cov <- matrix( Cov, d, d )
      dimnames(Cov) <- list( paste( "type", 1:d, sep="" ), paste( "type", 1:d, sep="" ) )
    }
    else
    {
      Cov <- NULL
      for( j in 1:d )
      {
        z0 <- rep(0,d)
	z0[j] <- 1
	Cov.j <- rowSums( V %*% diag( c( BGWM.mean( dists, type, d, (n-1), z0, maxiter ) ), d, d ) )
        for(i in (n-1):1)
        {
	  m.n_i <- m.n_i %*% m
	  if( i != 1 )
	    AUX <- rowSums( V %*% diag( c( BGWM.mean( dists, type, d, (i-1), z0, maxiter ) ), d, d ) )
	  else
	    AUX <- rowSums( V %*% diag( z0, d, d ) )
	  AUX <- matrix( AUX, d, d )
	  AUX <- t(m.n_i) %*% AUX %*% m.n_i
	  Cov.j <- Cov.j + AUX
        }
        Cov.j <- matrix( Cov.j, d, d )
	Cov <- rbind( Cov, Cov.j )
      }
      dimnames(Cov) <- list( paste( "dist", rep(1:d,rep(d,d)), ".type", rep(1:d,d), sep="" ), paste( "type", 1:d, sep="" ) )
    }
  }
  else
  {
    if( length(z0) != 0 )
    {
      V <- matrix( t(V), d*d, d )
      Cov <- rowSums( V %*% diag( z0, d, d ) )
      Cov <- matrix( Cov, d, d )
      dimnames(Cov) <- list( paste( "type", 1:d, sep="" ), paste( "type", 1:d, sep="" ) )
    }
    else
    {
      Cov <- V
      dimnames(Cov) <- list( paste( "dist", rep(1:d,rep(d,d)), ".type", rep(1:d,d), sep="" ), paste( "type", 1:d, sep="" ) )
    }
  }
  
  

  Cov

}



BGWM.gener.covar <- function(gener.dists, d)
{

  #Controles sobre los parametros

  s <- gener.dists[[1]]
  p <- gener.dists[[2]]
  v <- gener.dists[[3]]

  p <- data.frame( p = unlist(p), k = as.factor( rep( 1:d, unlist(s) ) ) )
  v <- data.frame( v = matrix( unlist( lapply( v, t ) ), ncol=d, byrow=TRUE ), k = as.factor( rep( 1:d, unlist(s) ) ) )
  #
  p[,1] <- p[,1] / rep( aggregate(p[,1], list(p[,2]) , sum)[,2] , unlist(s) )
  E2.X <- v
  E2.X[,-(d+1)] <- E2.X[,-(d+1)] * p[,1]
  E2.X <- aggregate( E2.X[,-(d+1)], list( E2.X[,(d+1)] ), sum )[,-1]
  E2.X <- matrix( c( apply( E2.X, 1, tcrossprod ) ), ncol=d, byrow=TRUE )
  # 
  aux1 <- v[,-(d+1)] * v[,-(d+1)] * p[,1]
  aux1 <- data.frame( aux1, k = as.factor( rep( 1:d, unlist(s) ) ) )
  aux1 <- aggregate( aux1[,-(d+1)], list(aux1[,(d+1)]), sum )[,-1]
  aux1 <- t( as.matrix( aux1 ) )
  # 
  a <- combn(1:d, 2)
  aux2 <- v[,a[1,]] * v[,a[2,]] * p[,1]
  aux2 <- data.frame( aux2, k = as.factor( rep( 1:d, unlist(s) ) ) )
  aux2 <- aggregate( aux2[,-(ncol(a)+1)], list(aux2[,(ncol(a)+1)]), sum )[,-1]
  aux2 <- t( as.matrix( aux2 ) )
  # 
  E.X2 <- matrix( rep( NA, (d*d*d) ), ncol=d*d )
  E.X2[row( E.X2 ) %% d == col( E.X2 ) %% d] <- aux1
  E.X2[( row( E.X2 ) %% d > col( E.X2 ) %% d | row( E.X2 ) == d ) & col( E.X2 ) %% d !=0] <- aux2
  E.X2[is.na(E.X2)] <- aux2[order(a[2,]),]

  covar <- t(E.X2) - E2.X

  covar

}



BGWM.multinom.covar <- function(multinom.dists, d, maxiter = 1e5)
{

  #Controles sobre los parametros

  dists <- multinom.dists[[1]]
  pmultinom <- multinom.dists[[2]]
  pmultinom <- as.matrix(pmultinom)
  mean <- rep( NA, d )
  var <- rep( NA, d )

  # unif
  a <- dists[,1] == "unif"
  if(TRUE %in% a)
  {
    mean[a] <- ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) ) / 2
    var[a] <- ( ( ( as.numeric( dists[a,3] ) - as.numeric( dists[a,2] ) + 1 ) ^ 2 ) - 1 ) / 12
  }
  # binom
  a <- dists[,1] == "binom"
  if(TRUE %in% a)
  {
    mean[a] <- as.numeric( dists[a,2] ) * as.numeric( dists[a,3] )
    var[a] <-  as.numeric( dists[a,2] ) * as.numeric( dists[a,3] ) * ( 1 - as.numeric( dists[a,3] ) )
  }
  # hyper
  a <- dists[,1] == "hyper"
  if(TRUE %in% a)
  {
    mean[a] <- as.numeric( dists[a,2] ) * as.numeric( dists[a,4] ) / ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) )
    var[a] <-  ( as.numeric( dists[a,2] ) * as.numeric( dists[a,4] ) * as.numeric( dists[a,3] ) * ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) - as.numeric( dists[a,4] ) ) ) / ( ( (as.numeric( dists[a,2] ) + as.numeric( dists[a,3] )) ^ 2 ) * ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) - 1 ) )
  }
  # geom
  a <- dists[,1] == "geom"
  if(TRUE %in% a)
  {
    mean[a] <- ( 1 - as.numeric( dists[a,2] ) ) / as.numeric( dists[a,2] )
    var[a] <- ( 1 -  as.numeric( dists[a,2] ) ) / ( as.numeric( dists[a,2] ) ^ 2 )
  }
  # nbinom
  a <- dists[,1] == "nbinom"
  if(TRUE %in% a)
  {
    mean[a] <- as.numeric( dists[a,2] ) * ( 1 - as.numeric( dists[a,3] ) ) / as.numeric( dists[a,3] )
    var[a] <-  as.numeric( dists[a,2] ) * ( 1 -  as.numeric( dists[a,3] ) ) / ( as.numeric( dists[a,3] ) ^ 2 )
  }
  # pois
  a <- dists[,1] == "pois"
  if(TRUE %in% a)
  {
    mean[a] <- as.numeric( dists[a,2] )
    var[a] <-  as.numeric( dists[a,2] )
  }
  # norm
  a <- dists[,1] == "norm"
  n <- length( var[a] )
  if(n > 0)
  {
    aux <- .C("param_estim_roundcut0_norm",
	       as.integer( maxiter ),
	       as.integer( n ),
	       as.double( as.numeric( dists[a,2] ) ),
	       as.double( as.numeric( dists[a,3] ) ),
	       as.integer( 3 ),
	       mean.estim=double( n ),
	       var.estim=double( n ),
	       PACKAGE="Branching")
    mean[a] <- round( aux$mean.estim, round( log10(maxiter)/2 ) )
    var[a] <- round( aux$var.estim, floor( log10(maxiter)/2 ) )
  }
  # lnorm
  a <- dists[,1] == "lnorm"
  n <- length( var[a] )
  if(n > 0)
  {
    aux <- .C("param_estim_round_lnorm",
	       as.integer( maxiter ),
	       as.integer( n ),
	       as.double( as.numeric( dists[a,2] ) ),
	       as.double( as.numeric( dists[a,3] ) ),
	       as.integer( 3 ),
	       mean.estim=double( n ),
	       var.estim=double( n ),
	       PACKAGE="Branching")
    mean[a] <- round( aux$mean.estim, round( log10(maxiter)/2 ) )
    var[a] <- round( aux$var.estim, floor( log10(maxiter)/2 ) )
  }
  # gamma
  a <- dists[,1] == "gamma"
  n <- length( var[a] )
  if(n > 0)
  {
    aux <- .C("param_estim_round_gamma",
	       as.integer( maxiter ),
	       as.integer( n ),
	       as.double( as.numeric( dists[a,2] ) ),
	       as.double( as.numeric( dists[a,3] ) ),
	       as.integer( 3 ),
	       mean.estim=double( n ),
	       var.estim=double( n ),
	       PACKAGE="Branching")
    mean[a] <- round( aux$mean.estim, round( log10(maxiter)/2 ) )
    var[a] <- round( aux$var.estim, floor( log10(maxiter)/2 ) )
  }

  aux1 <- matrix( rep( 0, (d*d*d) ), ncol=d )
  aux1[row(aux1) %% d == col(aux1) %% d] <- pmultinom
  aux2 <- matrix( c( apply( pmultinom, 1, tcrossprod ) ), ncol=d, byrow=TRUE )
  mean <- rep( mean, rep( d, d ) )
  var <- rep( var, rep( d, d ) )
  covar <- (aux1 - aux2) * mean + aux2 * var

  covar

}



BGWM.indep.covar <- function(indep.dists, d, maxiter = 1e5)
{

  #Controles sobre los parametros

  dists <- indep.dists
  var <- rep( NA, (d*d) )

  # unif
  a <- dists[,1] == "unif"
  if(TRUE %in% a)
    var[a] <- ( ( ( as.numeric( dists[a,3] ) - as.numeric( dists[a,2] ) + 1 ) ^ 2 ) - 1 ) / 12
  # binom
  a <- dists[,1] == "binom"
  if(TRUE %in% a)
    var[a] <-  as.numeric( dists[a,2] ) * as.numeric( dists[a,3] ) * ( 1 - as.numeric( dists[a,3] ) )
  # hyper
  a <- dists[,1] == "hyper"
  if(TRUE %in% a)
    var[a] <-  ( as.numeric( dists[a,2] ) * as.numeric( dists[a,4] ) * as.numeric( dists[a,3] ) * ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) - as.numeric( dists[a,4] ) ) ) / ( ( (as.numeric( dists[a,2] ) + as.numeric( dists[a,3] )) ^ 2 ) * ( as.numeric( dists[a,2] ) + as.numeric( dists[a,3] ) - 1 ) )
  # geom
  a <- dists[,1] == "geom"
  if(TRUE %in% a)
    var[a] <- ( 1 -  as.numeric( dists[a,2] ) ) / ( as.numeric( dists[a,2] ) ^ 2 )
  # nbinom
  a <- dists[,1] == "nbinom"
  if(TRUE %in% a)
    var[a] <-  as.numeric( dists[a,2] ) * ( 1 -  as.numeric( dists[a,3] ) ) / ( as.numeric( dists[a,3] ) ^ 2 )
  # pois
  a <- dists[,1] == "pois"
  if(TRUE %in% a)
    var[a] <-  as.numeric( dists[a,2] )
  # norm
  a <- dists[,1] == "norm"
  n <- length( var[a] )
  if(n > 0)
    var[a] <- round( .C("param_estim_roundcut0_norm",
			 as.integer( maxiter ),
			 as.integer( n ),
			 as.double( as.numeric( dists[a,2] ) ),
			 as.double( as.numeric( dists[a,3] ) ),
			 as.integer( 2 ),
			 mean.estim=double( n ),
			 var.estim=double( n ),
			 PACKAGE="Branching")$var.estim , floor( log10(maxiter)/2 ) )
  # lnorm
  a <- dists[,1] == "lnorm"
  n <- length( var[a] )
  if(n > 0)
    var[a] <- round( .C("param_estim_round_lnorm",
			 as.integer( maxiter ),
			 as.integer( n ),
			 as.double( as.numeric( dists[a,2] ) ),
			 as.double( as.numeric( dists[a,3] ) ),
			 as.integer( 2 ),
			 mean.estim=double( n ),
			 var.estim=double( n ),
			 PACKAGE="Branching")$var.estim , floor( log10(maxiter)/2 ) )
  # gamma
  a <- dists[,1] == "gamma"
  n <- length( var[a] )
  if(n > 0)
    var[a] <- round( .C("param_estim_round_gamma",
			 as.integer( maxiter ),
			 as.integer( n ),
			 as.double( as.numeric( dists[a,2] ) ),
			 as.double( as.numeric( dists[a,3] ) ),
			 as.integer( 2 ),
			 mean.estim=double( n ),
			 var.estim=double( n ),
			 PACKAGE="Branching")$var.estim , floor( log10(maxiter)/2 ) )

  covar <- matrix( rep( 0, (d*d*d) ), ncol=d )
  covar[row( covar ) %% d == col( covar ) %% d] <- matrix( var, ncol=d, byrow=TRUE )

  covar

}