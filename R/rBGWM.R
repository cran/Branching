# These functions simulate Bienayme - Galton - Watson multitype processes.
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


rBGWM <- function(dists, type=c("general","multinomial","independents"), d, n, z0=rep(1,d),
		  c.s=TRUE, tt.s=TRUE, rf.s=TRUE, file=NULL)
{

  type <- match.arg(type)
  type <- switch(type,
                 "general" = 1,
                 "multinomial" = 2,
                 "independents" = 3)

  R <- switch(type,
		 {#1
		   R <- rBGWM.gener(dists, d, n, z0, file)
		   R
		 },
		 {#2
		   R <- rBGWM.multinom(dists, d, n, z0, file)
		   R
		 },
		 {#3
		   R <- rBGWM.indep(dists, d, n, z0, file)
		   R
		 })

  cdata <- R$cdata
  if(TRUE %in% (cdata<0))
    warning("exceeded maximum capacity of data type. Process truncated")
  cdata <- matrix(as.numeric(cdata),ncol=d,byrow=TRUE)
  cdata <- data.frame(cdata)
  dimnames(cdata) <- list(paste("i",
                                rep(1:n,rep(d,n)),
                                ".",
                                rep(paste("type",1:d,sep=""),n),sep=""),
                          paste("type",1:d,sep=""))
  cdata$i <- as.factor(rep(1:n,rep(d,n)))

  ttdata <- aggregate(cdata[,1:d],list(cdata$i),sum)[,-1]
  ttdata <- as.matrix(ttdata)
  cdata <- as.matrix(cdata[,1:d])

  if(rf.s==TRUE)
  {
    rfdata <- ttdata / apply( ttdata, 1, sum )
    rfdata <- as.matrix(rfdata)
  }
  else
    rfdata <- NULL

  if(tt.s==FALSE)
    ttdata <- NULL

  if(c.s==FALSE)
    cdata <- NULL

  out <- list(i.dists=dists,
              i.d=d,
              i.n=n,
              i.z0=z0,
              o.c.s=cdata,
              o.tt.s=ttdata,
              o.rf.s=rfdata)

  out

}



# rBGWM general
rBGWM.gener <- function(gener.dists, d, n, z0=rep(1,d), file=NULL)
{

  sizes <- c( unlist( gener.dists[[1]] ) )
  probs <- lapply( gener.dists[[2]], cumsum )
  probs <- c( unlist( probs ) )
  vectors <- lapply( gener.dists[[3]], unique )
  vectors <- c( unlist( lapply( vectors, t ) ) )

  if(length(d)!= 1)
    stop("'d' must be a positive number")
  if(length(n)!= 1)
    stop("'n' must be a positive number")
  if(length(z0)!= d)
    stop("'z0' must be a d-dimensional vector")
  if(TRUE %in% (z0 < 0))
    stop("'z0' must have positive elements")
  if(length(sizes) != d)
    stop("'gener.dists$sizes' must be a d-dimensional vector")
  if(length(probs) != sum(sizes))
    stop("'gener.dists$probs' does not have the right structure (wrong number of elements)")
  if(TRUE %in% (probs <= 0))
    stop("'gener.dists$probs' elements must be all positive")
  if(length(vectors) != sum(sizes*d))
    stop("'gener.dists$vectors' does not have the right structure (wrong number of elements or duplicated rows by distribution)")
  # more restrictions?

  R <- .C(rBGWMgeneral,
          as.integer(d),
          as.integer(n),
          as.integer(z0),
          as.integer(sizes),
          as.integer(vectors),
          as.double(probs),
          cdata=double(d*d*n),
	  as.character(file))

  R

}



# rBGWM multinomial
rBGWM.multinom <- function(multinom.dists, d, n, z0=rep(1,d), file=NULL)
{

  dists <- multinom.dists[[1]]
  pmultinom <- multinom.dists[[2]]
  p <- as.matrix(pmultinom)
  p <- p / apply( p, 1, sum )
  nrodists <- nrow(dists)
  dists[is.na(dists)] <- 0
  names.dists <- dists[,1]
  param.dists <- as.matrix(dists[,-1])
  aux <- d*d*n

  if(nrow(p) != ncol(p) || nrow(p) != d || ncol(p) != d)
    stop("'pmultinom' must be a squared matrix of order d")
  if(length(d) != 1)
    stop("'d' must be a number")
  if(length(n) != 1)
    stop("'n' must be a number")
  if(length(z0) != d)
    stop("'z0' must be a d-dimensional vector")
  if(TRUE %in% (z0 < 0))
    stop("'z0' must have positive elements")
  if(nrodists != 1 && nrodists != d)
    stop("'dists' must be 1 or d distributions")
  if(FALSE %in% (tolower(names.dists) %in% c("binom","gamma","geom","hyper","lnorm","nbinom","norm","pois","unif")))
    stop("There are in 'dists' some distributions not implemented yet (only binom, gamma, geom, hyper, lnorm, nbinom, norm, pois, unif are implemented)")
  # more restrictions?

  names.dists[names.dists == "unif"] <- 1
  names.dists[names.dists == "binom"] <- 2
  names.dists[names.dists == "hyper"] <- 3
  names.dists[names.dists == "geom"] <- 4
  names.dists[names.dists == "nbinom"] <- 5
  names.dists[names.dists == "pois"] <- 6
  names.dists[names.dists == "norm"] <- 7
  names.dists[names.dists == "lnorm"] <- 8
  names.dists[names.dists == "gamma"] <- 9

  R <- .C(rBGWMmultinomial,
          as.integer(d),
          as.integer(n),
          as.integer(z0),
          as.integer(nrodists),
          as.integer(names.dists),
          as.integer(ncol(param.dists)),
          as.double(t(param.dists)),
          as.double(t(p)),
          cdata=double(aux),
	  as.character(file))

  R

}



# rBGWM independents
rBGWM.indep <- function(indep.dists, d, n, z0=rep(1,d), file=NULL)
{

  dists <- indep.dists
  nrodists <- nrow(dists)
  names.dists <- dists[,1]
  dists[is.na(dists)] <- 0
  param.dists <- as.matrix(dists[,-1])
  aux <- d*d*n

  if(length(d)!= 1)
    stop("'d' must be a number")
  if(length(n)!= 1)
    stop("'n' must be a number")
  if(length(z0)!= d)
    stop("'z0' must be a d-dimensional vector")
  if(TRUE %in% (z0 < 0))
    stop("'z0' must have positive elements")
  if(nrodists != 1 && nrodists != (d*d))
    stop("'dists' must be 1 or d*d distributions")
  if(FALSE %in% (tolower(names.dists) %in% c("binom","gamma","geom","hyper","lnorm","nbinom","norm","pois","unif")))
    stop("There are in 'dists' some distributions not implemented yet (only binom, gamma, geom, hyper, lnorm, nbinom, norm, pois, unif are implemented)")
  # more restrictions?

  names.dists[names.dists == "unif"] <- 1
  names.dists[names.dists == "binom"] <- 2
  names.dists[names.dists == "hyper"] <- 3
  names.dists[names.dists == "geom"] <- 4
  names.dists[names.dists == "nbinom"] <- 5
  names.dists[names.dists == "pois"] <- 6
  names.dists[names.dists == "norm"] <- 7
  names.dists[names.dists == "lnorm"] <- 8
  names.dists[names.dists == "gamma"] <- 9

  R <- .C(rBGWMindependent,
          as.integer(d),
          as.integer(n),
          as.integer(z0),
          as.integer(nrodists),
          as.integer(names.dists),
          as.integer(ncol(param.dists)),
          as.double(t(param.dists)),
          cdata=double(aux),
	  as.character(file))

  R

}
