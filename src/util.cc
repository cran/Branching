/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

//

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <cstdio>
#include "util.h"

using namespace std;

extern "C"  {

  void param_estim_roundcut0_norm( int *max_iter, int *n_dists, double *param1, double *param2, 
				   int *n_estim ,double *mean_estim, double *var_estim )
  {
    int max = *max_iter;
    int size = *n_dists;
    double mean=0.0, var=0.0;
    int n = *n_estim;
    int i, l;
    int aux;

    GetRNGstate();

    switch(n)
    {
      case 1:
	for( i = 0 ; i < size ; i++ )
	{
	  for( l = 1 ; l <= max ; l++ )
          {
            mean += cut0( fround( rnorm( param1[i], param2[i] ) , 0) );
          } // end for(l)
	  mean /= max;
          mean_estim[i] = mean;
        } // end for(i)
        break;
      case 2:
	for( i = 0 ; i < size ; i++ )
	{
          for( l = 1 ; l <= max ; l++ )
          {
            aux = cut0( fround( rnorm( param1[i], param2[i] ) , 0) );
	    mean += aux;
	    var += aux * aux;
          } // end for(l)
	  mean /= max;
	  var /= max;
	  var -= (mean) * (mean);
          var_estim[i] = var;
	} // end for(i)
	break;
      case 3:
	for( i = 0 ; i < size ; i++ )
	{
          for( l = 1 ; l <= max ; l++ )
          {
            aux = cut0( fround( rnorm( param1[i], param2[i] ) , 0) );
	    mean += aux;
	    var += aux * aux;
          } // end for(l)
	  mean /= max;
	  var /= max;
	  var -= (mean) * (mean);
          var_estim[i] = var;
	  mean_estim[i] = mean;
	} // end for(i)
	break;
    } //end switch

    PutRNGstate();

  }
  
  void param_estim_round_lnorm( int *max_iter, int *n_dists, double *param1, double *param2, 
				int *n_estim ,double *mean_estim, double *var_estim )
  {
    int max = *max_iter;
    int size = *n_dists;
    double mean=0.0, var=0.0;
    int n = *n_estim;
    int i, l;
    int aux;

    GetRNGstate();

    switch(n)
    {
      case 1:
	for( i = 0 ; i < size ; i++ )
	{
	  for( l = 1 ; l <= max ; l++ )
          {
            mean += fround( rlnorm( param1[i], param2[i] ) , 0);
          } // end for(l)
	  mean /= max;
          mean_estim[i] = mean;
        } // end for(i)
        break;
      case 2:
	for( i = 0 ; i < size ; i++ )
	{
          for( l = 1 ; l <= max ; l++ )
          {
            aux = fround( rlnorm( param1[i], param2[i] ) , 0);
	    mean += aux;
	    var += aux * aux;
          } // end for(l)
	  mean /= max;
	  var /= max;
	  var -= (mean) * (mean);
          var_estim[i] = var;
	} // end for(i)
	break;
      case 3:
	for( i = 0 ; i < size ; i++ )
	{
          for( l = 1 ; l <= max ; l++ )
          {
            aux = fround( rlnorm( param1[i], param2[i] ) , 0);
	    mean += aux;
	    var += aux * aux;
          } // end for(l)
	  mean /= max;
	  var /= max;
	  var -= (mean) * (mean);
          var_estim[i] = var;
	  mean_estim[i] = mean;
	} // end for(i)
	break;
    } //end switch

    PutRNGstate();

  }
  
  void param_estim_round_gamma( int *max_iter, int *n_dists, double *param1, double *param2, 
				int *n_estim ,double *mean_estim, double *var_estim )
  {
    int max = *max_iter;
    int size = *n_dists;
    double mean=0.0, var=0.0;
    int n = *n_estim;
    int i, l;
    int aux;

    GetRNGstate();

    switch(n)
    {
      case 1:
	for( i = 0 ; i < size ; i++ )
	{
	  for( l = 1 ; l <= max ; l++ )
          {
            mean += fround( rgamma( param1[i], param2[i] ) , 0);
          } // end for(l)
	  mean /= max;
          mean_estim[i] = mean;
        } // end for(i)
        break;
      case 2:
	for( i = 0 ; i < size ; i++ )
	{
          for( l = 1 ; l <= max ; l++ )
          {
            aux = fround( rgamma( param1[i], param2[i] ) , 0);
	    mean += aux;
	    var += aux * aux;
          } // end for(l)
	  mean /= max;
	  var /= max;
	  var -= (mean) * (mean);
          var_estim[i] = var;
	} // end for(i)
	break;
      case 3:
	for( i = 0 ; i < size ; i++ )
	{
          for( l = 1 ; l <= max ; l++ )
          {
            aux = fround( rgamma( param1[i], param2[i] ) , 0);
	    mean += aux;
	    var += aux * aux;
          } // end for(l)
	  mean /= max;
	  var /= max;
	  var -= (mean) * (mean);
          var_estim[i] = var;
	  mean_estim[i] = mean;
	} // end for(i)
	break;
    } //end switch

    PutRNGstate();
    
  }

} //end extern
