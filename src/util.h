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

#ifndef UTIL_H
#define UTIL_H
extern "C" 
{
  void param_estim_roundcut0_norm( int *max_iter, int *n_dists, double *param1, double *param2, 
				   int *n_estim ,double *mean_estim, double *var_estim );
  void param_estim_round_lnorm( int *max_iter, int *n_dists, double *param1, double *param2, 
				int *n_estim ,double *mean_estim, double *var_estim );
  void param_estim_round_gamma( int *max_iter, int *n_dists, double *param1, double *param2, 
				int *n_estim ,double *mean_estim, double *var_estim );
}
  inline int cut0(int a){return a < 0 ? 0 : a;};
  inline double cut0(double a){return a < 0 ? 0 : a;};
#endif
