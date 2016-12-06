/*
# These functions simulate Bienayme - Galton - Watson multitype processes. 
# Copyright (C) 2010  Camilo Jos? Torres Jim?nez <cjtorresj@unal.edu.co>

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
*/

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <cstdio>
#include "rBGWM.h"
#include "util.h"

using namespace std;

extern "C"  {

// rBGWMmultinomial function
  void rBGWMmultinomial(int *d, int *n, unsigned long *z0, int *nrodists, int *names_dists, int *nroparam_dists, 
                        double *param_dists, double *pmultinom, double *cdata, char **outfile)
  {
    
    int nro_types = *d; 
    int max_time = *n; 
    int nro_dists = *nrodists; 
    int nro_param = *nroparam_dists;
    int i, j, k, m = 0, o;
    unsigned long l;
    int aux, aux2, nro_hijos;
    int *result;
    FILE *of = NULL;

    if(outfile != NULL)
    {
      of = fopen(*outfile, "w");
      if(of == NULL)
      {
	//printf("Failure to open the file\n");
	return;
      }
      for(o = 0 ; o < nro_types ; o++)
	fprintf(of,"\ttype%d",o+1);
      fprintf(of,"\n");
      fflush(of);
    }
    
    if((result = new int[nro_types]) == NULL)
    {
      //printf("Out of memory\n");
      return;
    }

    GetRNGstate();
    
    //
    for(k = 0 ; k < nro_types ; k++)
    {
      //
      for(o = 0 ; o < nro_types ; o++)
        cdata[o + nro_types * k] = 0;
      //
      if(z0[k] > 0)
      {
        //
        if(nro_dists != 1)
          m = k;
        aux2 = nro_param * m;
        //
        switch(names_dists[m])
        {
          //unif=1,binom=2,hyper=3,geom=4,nbinom=5,pois=6,norm=7,lnorm=8,gamma=9
          case 1:
            //
            for(l = 1 ; l <= z0[k] ; l++)
            {
              //
              nro_hijos = (int) runif(param_dists[aux2], param_dists[1 + aux2] + 1);
              //
              if(nro_hijos > 0)
              {
                //
                rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                //
                for(o = 0 ; o < nro_types ; o++)
                  cdata[o + nro_types * k] += result[o];
              } //end if
            } //end for(l)
            break;
          case 2:
            //
            for(l = 1 ; l <= z0[k] ; l++)
            {
              //
              nro_hijos = rbinom(param_dists[aux2], param_dists[1 + aux2]);
              //
              if(nro_hijos > 0)
              {
                //
                rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                //
                for(o = 0 ; o < nro_types ; o++)
                  cdata[o + nro_types * k] += result[o];
              } //end if
            } //end for(l)
            break;
          case 3:
            //
            for(l = 1 ; l <= z0[k] ; l++)
            {
              //
              nro_hijos = rhyper(param_dists[aux2], param_dists[1 + aux2], param_dists[2 + aux2]);
              //
              if(nro_hijos > 0)
              {
                //
                rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                //
                for(o = 0 ; o < nro_types ; o++)
                  cdata[o + nro_types * k] += result[o];
              } //end if
            } //end for(l)
            break;
          case 4:
            //
            for(l = 1 ; l <= z0[k] ; l++)
            {
              //
              nro_hijos = rgeom(param_dists[aux2]);
              //
              if(nro_hijos > 0)
              {
                //
                rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                //
                for(o = 0 ; o < nro_types ; o++)
                  cdata[o + nro_types * k] += result[o];
              } //end if
            } //end for(l)
            break;
          case 5:
            //
            for(l = 1 ; l <= z0[k] ; l++)
            {
              //
              nro_hijos = rnbinom(param_dists[aux2], param_dists[1 + aux2]);
              //
              if(nro_hijos > 0)
              {
                //
                rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                //
                for(o = 0 ; o < nro_types ; o++)
                  cdata[o + nro_types * k] += result[o];
              } //end if
            } //end for(l)
            break;
          case 6:
            //
            for(l = 1 ; l <= z0[k] ; l++)
            {
              //
              nro_hijos = rpois(param_dists[aux2]);
              //
              if(nro_hijos > 0)
              {
                //
                rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                //
                for(o = 0 ; o < nro_types ; o++)
                  cdata[o + nro_types * k] += result[o];
              } //end if
            } //end for(l)
            break;
          case 7:
            //
            for(l = 1 ; l <= z0[k] ; l++)
            {
              //
              nro_hijos = cut0( fround(rnorm(param_dists[aux2], param_dists[1 + aux2]),0) );
              //
              if(nro_hijos > 0)
              {
                //
                rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                //
                for(o = 0 ; o < nro_types ; o++)
                  cdata[o + nro_types * k] += result[o];
              } //end if
            } //end for(l)
            break;
          case 8:
            //
            for(l = 1 ; l <= z0[k] ; l++)
            {
              //
              nro_hijos = fround(rlnorm(param_dists[aux2], param_dists[1 + aux2]),0);
              //
              if(nro_hijos > 0)
              {
                //
                rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                //
                for(o = 0 ; o < nro_types ; o++)
                  cdata[o + nro_types * k] += result[o];
              } //end if
            } //end for(l)
            break;
          case 9:
            //
            for(l = 1 ; l <= z0[k] ; l++)
            {
              //
              nro_hijos = fround(rgamma(param_dists[aux2], param_dists[1 + aux2]),0);
              //
              if(nro_hijos > 0)
              {
                //
                rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                //
                for(o = 0 ; o < nro_types ; o++)
                  cdata[o + nro_types * k] += result[o];
              } //end if
            } //end for(l)
            break;
        } //end switch
      } //end if
      else
       if(z0[k] < 0)
         return;
      if(outfile != NULL)
      {
	fprintf(of,"i1.type%d\t",k+1);
	for(o = 0 ; o < nro_types ; o++)
	   fprintf(of,"%.0f\t",cdata[o + nro_types * k]);
	fprintf(of,"\n");
	fflush(of);
      }
    } //end for(k)

    //
    for(i = 1 ; i < max_time ; i++)
    {
      //
      for(k = 0 ; k < nro_types ; k++)
      {
        for(o = 0 ; o < nro_types ; o++)
          cdata[o + nro_types * k + nro_types * nro_types * i] = 0;
        //
        if(nro_dists != 1)
          m = k;
        aux2 = nro_param * m;
        //
        switch(names_dists[m])
        {
          //unif=1,binom=2,hyper=3,geom=4,nbinom=5,pois=6,norm=7,lnorm=8,gamma=9
          case 1:
            //
            for(j = 0 ; j < nro_types ; j++)
            {
              //
              aux = k + nro_types * j + nro_types * nro_types * (i-1);
              if(cdata[aux] > 0)
              {
                //
                for(l = 1 ; l <= cdata[aux] ; l++)
                {
                  nro_hijos = (int) runif(param_dists[aux2], param_dists[1 + aux2] + 1);
                  //
                  if(nro_hijos > 0)
                  {
                    //
                    rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                    //
                    for(o = 0 ; o < nro_types ; o++)
                      cdata[o + nro_types * k + nro_types * nro_types * i] += result[o];
                  } //if
                } //end for(l)
              } //end if
              else
                if(cdata[aux] < 0)
                  return;
            } //end for(j)
            break;
          case 2:
            //
            for(j = 0 ; j < nro_types ; j++)
            {
              //
              aux = k + nro_types * j + nro_types * nro_types * (i-1);
              if(cdata[aux] > 0)
              {
                //
                for(l = 1 ; l <= cdata[aux] ; l++)
                {
                  nro_hijos = rbinom(param_dists[aux2], param_dists[1 + aux2]);
                  //
                  if(nro_hijos > 0)
                  {
                    //
                    rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                    //
                    for(o = 0 ; o < nro_types ; o++)
                      cdata[o + nro_types * k + nro_types * nro_types * i] += result[o];
                  } //if
                } //end for(l)
              } //end if
              else
                if(cdata[aux] < 0)
                  return;
            } //end for(j)
            break;
          case 3:
            //
            for(j = 0 ; j < nro_types ; j++)
            {
              //
              aux = k + nro_types * j + nro_types * nro_types * (i-1);
              if(cdata[aux] > 0)
              {
                //
                for(l = 1 ; l <= cdata[aux] ; l++)
                {
                  nro_hijos = rhyper(param_dists[aux2], param_dists[1 + aux2], param_dists[2 + aux2]);
                  //
                  if(nro_hijos > 0)
                  {
                    //
                    rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                    //
                    for(o = 0 ; o < nro_types ; o++)
                      cdata[o + nro_types * k + nro_types * nro_types * i] += result[o];
                  } //if
                } //end for(l)
              } //end if
              else
                if(cdata[aux] < 0)
                  return;
            } //end for(j)
            break;
          case 4:
            //
            for(j = 0 ; j < nro_types ; j++)
            {
              //
              aux = k + nro_types * j + nro_types * nro_types * (i-1);
              if(cdata[aux] > 0)
              {
                //
                for(l = 1 ; l <= cdata[aux] ; l++)
                {
                  nro_hijos = rgeom(param_dists[aux2]);
                  //
                  if(nro_hijos > 0)
                  {
                    //
                    rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                    //
                    for(o = 0 ; o < nro_types ; o++)
                      cdata[o + nro_types * k + nro_types * nro_types * i] += result[o];
                  } //if
                } //end for(l)
              } //end if
              else
                if(cdata[aux] < 0)
                  return;
            } //end for(j)
            break;
          case 5:
            //
            for(j = 0 ; j < nro_types ; j++)
            {
              //
              aux = k + nro_types * j + nro_types * nro_types * (i-1);
              if(cdata[aux] > 0)
              {
                //
                for(l = 1 ; l <= cdata[aux] ; l++)
                {
                  nro_hijos = rnbinom(param_dists[aux2], param_dists[1 + aux2]);
                  //
                  if(nro_hijos > 0)
                  {
                    //
                    rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                    //
                    for(o = 0 ; o < nro_types ; o++)
                      cdata[o + nro_types * k + nro_types * nro_types * i] += result[o];
                  } //if
                } //end for(l)
              } //end if
              else
                if(cdata[aux] < 0)
                  return;
            } //end for(j)
            break;
          case 6:
            //
            for(j = 0 ; j < nro_types ; j++)
            {
              //
              aux = k + nro_types * j + nro_types * nro_types * (i-1);
              if(cdata[aux] > 0)
              {
                //
                for(l = 1 ; l <= cdata[aux] ; l++)
                {
                  nro_hijos = rpois(param_dists[aux2]);
                  //
                  if(nro_hijos > 0)
                  {
                    //
                    rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                    //
                    for(o = 0 ; o < nro_types ; o++)
                      cdata[o + nro_types * k + nro_types * nro_types * i] += result[o];
                  } //if
                } //end for(l)
              } //end if
              else
                if(cdata[aux] < 0)
                  return;
            } //end for(j)
            break;
          case 7:
            //
            for(j = 0 ; j < nro_types ; j++)
            {
              //
              aux = k + nro_types * j + nro_types * nro_types * (i-1);
              if(cdata[aux] > 0)
              {
                //
                for(l = 1 ; l <= cdata[aux] ; l++)
                {
                  nro_hijos = cut0( fround(rnorm(param_dists[aux2], param_dists[1 + aux2]),0) );
                  //
                  if(nro_hijos > 0)
                  {
                    //
                    rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                    //
                    for(o = 0 ; o < nro_types ; o++)
                      cdata[o + nro_types * k + nro_types * nro_types * i] += result[o];
                  } //if
                } //end for(l)
              } //end if
              else
                if(cdata[aux] < 0)
                  return;
            } //end for(j)
            break;
          case 8:
            //
            for(j = 0 ; j < nro_types ; j++)
            {
              //
              aux = k + nro_types * j + nro_types * nro_types * (i-1);
              if(cdata[aux] > 0)
              {
                //
                for(l = 1 ; l <= cdata[aux] ; l++)
                {
                  nro_hijos = fround(rlnorm(param_dists[aux2], param_dists[1 + aux2]),0);
                  //
                  if(nro_hijos > 0)
                  {
                    //
                    rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                    //
                    for(o = 0 ; o < nro_types ; o++)
                      cdata[o + nro_types * k + nro_types * nro_types * i] += result[o];
                  } //if
                } //end for(l)
              } //end if
              else
                if(cdata[aux] < 0)
                  return;
            } //end for(j)
            break;
          case 9:
            //
            for(j = 0 ; j < nro_types ; j++)
            {
              //
              aux = k + nro_types * j + nro_types * nro_types * (i-1);
              if(cdata[aux] > 0)
              {
                //
                for(l = 1 ; l <= cdata[aux] ; l++)
                {
                  nro_hijos = fround(rgamma(param_dists[aux2], param_dists[1 + aux2]),0);
                  //
                  if(nro_hijos > 0)
                  {
                    //
                    rmultinom(nro_hijos, &pmultinom[nro_types * k], nro_types, result);
                    //
                    for(o = 0 ; o < nro_types ; o++)
                      cdata[o + nro_types * k + nro_types * nro_types * i] += result[o];
                  } //if
                } //end for(l)
              } //end if
              else
                if(cdata[aux] < 0)
                  return;
            } //end for(j)
            break;
        } //end switch
	if(outfile != NULL)
	{
	  fprintf(of, "i%d.type%d\t", i+1, k+1);
	  for(o = 0 ; o < nro_types ; o++)
	     fprintf(of, "%.0f\t", cdata[o + nro_types * k + nro_types * nro_types * i]);
	  fprintf(of,"\n");
	  fflush(of);
	}
      } //end for(k)
    } //end for(i)

    PutRNGstate();
    
    if(outfile != NULL)
    {
      fclose(of);
    }
    
  } //end rBGWMmultinom


// rBGWMgeneral function
  void rBGWMgeneral(int *d, int *n, unsigned long *z0, int *sizes, int *vectors, double *probs, 
                    double *cdata, char **outfile)
  {

    int nro_types = *d; 
    int max_time = *n; 
    int i, j, k, o, m;
    unsigned long l;
    int aux_size = 0, aux_size_n, aux_i_1;
    double aux_runif;
    FILE *of = NULL;

    if(outfile != NULL)
    {
      of = fopen(*outfile, "w");
      if(of == NULL)
      {
	//printf("Failure to open the file\n");
	return;
      }
      for(o = 0 ; o < nro_types ; o++)
	fprintf(of,"\ttype%d",o+1);
      fprintf(of,"\n");
      fflush(of);
    }
    
    GetRNGstate();
    
    //
    for(k = 0 ; k < nro_types ; k++)
    {
      //
      aux_size_n = aux_size + sizes[k];
      for(o = 0 ; o < aux_size_n ; o++)
        probs[o + aux_size] /= probs[aux_size_n - 1];
      //
      for(o = 0 ; o < nro_types ; o++)
        cdata[o + k * nro_types] = 0;
      //
      if(z0[k] > 0)
      {
        //
        for(l = 0 ; l < z0[k] ; l++)
        {
          //
          aux_runif = runif(0, 1);
          //
          for(m = 0 ; probs[m + aux_size] < aux_runif ; m++){}
          //
          for(o = 0 ; o < nro_types ; o++)
          {
            cdata[o + k * nro_types] += vectors[o + m * nro_types + aux_size * nro_types];
          } //end for(o)
        } //end for(l)
      } //end if
      else
        if(z0[k] < 0)
          return;
      //
      aux_size = aux_size_n;
      if(outfile != NULL)
      {
	fprintf(of,"i1.type%d\t",k+1);
	for(o = 0 ; o < nro_types ; o++)
	   fprintf(of,"%.0f\t",cdata[o + nro_types * k]);
	fprintf(of,"\n");
	fflush(of);
      }
    } //end for(k)

    //
    for(i = 1 ; i < max_time ; i++)
    {
      //
      aux_size = 0;
      //
      for(k = 0 ; k < nro_types ; k++)
      {
        //
        for(o = 0 ; o < nro_types ; o++)
          cdata[o + k * nro_types + i * nro_types * nro_types] = 0;
        //
        for(j = 0 ; j < nro_types ; j++)
        {
          //
          aux_i_1 = k + j * nro_types + (i - 1) * nro_types * nro_types;
          //
          if(cdata[aux_i_1] > 0)
          {
            //
            for(l = 0 ; l < cdata[aux_i_1] ; l++)
            {
              //
              aux_runif = runif(0, 1);
              //
              for(m = 0 ; probs[m + aux_size] < aux_runif ; m++){}
              //
              for(o = 0 ; o < nro_types ; o++)
              {
                cdata[o + k * nro_types + i * nro_types * nro_types] += vectors[o + m * nro_types + aux_size * nro_types];
              } //end for(o)
            } //end for(l)
          } //end if
          else
            if(cdata[aux_i_1] < 0)
              return;
        } //end for(j)
        //
        aux_size = aux_size + sizes[k];
	if(outfile != NULL)
	{
	  fprintf(of, "i%d.type%d\t", i+1, k+1);
	  for(o = 0 ; o < nro_types ; o++)
	     fprintf(of, "%.0f\t", cdata[o + nro_types * k + nro_types * nro_types * i]);
	  fprintf(of,"\n");
	  fflush(of);
	}
      } //end for(k)
    } //end for(i)    

    PutRNGstate();

    if(outfile != NULL)
    {
      fclose(of);
    }
      
  } //end rBGWMgeneral


// rBGWMindependent function
  void rBGWMindependent(int *d, int *n, unsigned long *z0, int *nrodists, int *names_dists, int *nroparam_dists, 
                        double *param_dists, double *cdata, char **outfile)
  {
    
    int nro_types = *d; 
    int max_time = *n; 
    int nro_dists = *nrodists; 
    int nro_param = *nroparam_dists;
    int i, j, k, m = 0, o;
    unsigned long l;
    int aux, aux2, aux_i, aux_i_1;
    FILE *of = NULL;

    if(outfile != NULL)
    {
      of = fopen(*outfile, "w");
      if(of == NULL)
      {
	//printf("Failure to open the file\n");
	return;
      }
      for(o = 0 ; o < nro_types ; o++)
	fprintf(of,"\ttype%d",o+1);
      fprintf(of,"\n");
      fflush(of);
    }
    
    GetRNGstate();
    
    //
    for(k = 0 ; k < nro_types ; k++)
    {
      //
      if(z0[k] > 0)
      {
        //
        if(nro_dists != 1)
          m = k;
        //
        for(o = 0 ; o < nro_types ; o++)
        {
          //
          aux = o + nro_types * k;
          cdata[aux] = 0;
          aux2 = nro_param * (o + nro_types * m);
          //
          switch(names_dists[o + nro_types * m])
          {
            //unif=1,binom=2,hyper=3,geom=4,nbinom=5,pois=6,norm=7,lnorm=8,gamma=9
            case 1:
              //
              for(l = 1 ; l <= z0[k] ; l++)
              {
                //
                cdata[aux] += (int) runif(param_dists[aux2], param_dists[1 + aux2] + 1);
              } //end for(l)
              break;
            case 2:
              //
              for(l = 1 ; l <= z0[k] ; l++)
              {
                //
                cdata[aux] += rbinom(param_dists[aux2], param_dists[1 + aux2]);
              } //end for(l)
              break;
            case 3:
              //
              for(l = 1 ; l <= z0[k] ; l++)
              {
                //
                cdata[aux] += rhyper(param_dists[aux2], param_dists[1 + aux2], param_dists[2 + aux2]);
              } //end for(l)
              break;
            case 4:
              //
              for(l = 1 ; l <= z0[k] ; l++)
              {
                //
                cdata[aux] += rgeom(param_dists[aux2]);
              } //end for(l)
              break;
            case 5:
              //
              for(l = 1 ; l <= z0[k] ; l++)
              {
                //
                cdata[aux] += rnbinom(param_dists[aux2], param_dists[1 + aux2]);
              } //end for(l)
              break;
            case 6:
              //
              for(l = 1 ; l <= z0[k] ; l++)
              {
                //
                cdata[aux] += rpois(param_dists[aux2]);
              } //end for(l)
              break;
            case 7:
              //
              for(l = 1 ; l <= z0[k] ; l++)
              {
                //
                cdata[aux] += cut0( fround(rnorm(param_dists[aux2], param_dists[1 + aux2]),0) );
              } //end for(l)
              break;
            case 8:
              //
              for(l = 1 ; l <= z0[k] ; l++)
              {
                //
                cdata[aux] += fround(rlnorm(param_dists[aux2], param_dists[1 + aux2]),0);
              } //end for(l)
              break;
            case 9:
              //
              for(l = 1 ; l <= z0[k] ; l++)
              {
                //
                cdata[aux] += fround(rgamma(param_dists[aux2], param_dists[1 + aux2]),0);
              } //end for(l)
              break;
          } //end switch
        } //end for(o)
      } //end if
      else
        if(z0[k] < 0)
          return;
      if(outfile != NULL)
      {
	fprintf(of,"i1.type%d\t",k+1);
	for(o = 0 ; o < nro_types ; o++)
	   fprintf(of,"%.0f\t",cdata[o + nro_types * k]);
	fprintf(of,"\n");
	fflush(of);
      }
    } //end for(k)

    //
    for(i = 1 ; i < max_time ; i++)
    {
      //
      for(k = 0 ; k < nro_types ; k++)
      {
        //
        if(nro_dists != 1)
          m = k;
        //
        for(o = 0 ; o < nro_types ; o++)
        {
          aux = o + nro_types * m;
          aux_i = o + nro_types * k + nro_types * nro_types * i;
          cdata[aux_i] = 0;
          //
          for(j = 0 ; j < nro_types ; j++)
          {
            //
            aux_i_1 = k + nro_types * j + nro_types * nro_types * (i-1);
            aux2 = nro_param * aux;
            //
            if(cdata[aux_i_1] > 0)
            {
              //
              switch(names_dists[aux])
              {
                //unif=1,binom=2,hyper=3,geom=4,nbinom=5,pois=6,norm=7,lnorm=8,gamma=9
                case 1:
                  //
                  for(l = 1 ; l <= cdata[aux_i_1] ; l++)
                  {
                    cdata[aux_i] += (int) runif(param_dists[aux2], param_dists[1 + aux2] + 1);
                  } //end for(l)
                  break;
                case 2:
                  //
                  for(l = 1 ; l <= cdata[aux_i_1] ; l++)
                  {
                    cdata[aux_i] += rbinom(param_dists[aux2], param_dists[1 + aux2]);
                  } //end for(l)
                  break;
                case 3:
                  //
                  for(l = 1 ; l <= cdata[aux_i_1] ; l++)
                  {
                    cdata[aux_i] += rhyper(param_dists[aux2], param_dists[1 + aux2], param_dists[2 + aux2]);
                  } //end for(l)
                  break;
                case 4:
                  //
                  for(l = 1 ; l <= cdata[aux_i_1] ; l++)
                  {
                    cdata[aux_i] += rgeom(param_dists[aux2]);
                  } //end for(l)
                  break;
                case 5:
                  //
                  for(l = 1 ; l <= cdata[aux_i_1] ; l++)
                  {
                    cdata[aux_i] += rnbinom(param_dists[aux2], param_dists[1 + aux2]);
                  } //end for(l)
                  break;
                case 6:
                  //
                  for(l = 1 ; l <= cdata[aux_i_1] ; l++)
                  {
                    cdata[aux_i] += rpois(param_dists[aux2]);
                  } //end for(l)
                  break;
                case 7:
                  //
                  for(l = 1 ; l <= cdata[aux_i_1] ; l++)
                  {
                    cdata[aux_i] += cut0( fround(rnorm(param_dists[aux2], param_dists[1 + aux2]),0) );
                  } //end for(l)
                  break;
                case 8:
                  //
                  for(l = 1 ; l <= cdata[aux_i_1] ; l++)
                  {
                    cdata[aux_i] += fround(rlnorm(param_dists[aux2], param_dists[1 + aux2]),0);
                  } //end for(l)
                  break;
               case 9:
                  //
                  for(l = 1 ; l <= cdata[aux_i_1] ; l++)
                  {
                    cdata[aux_i] += fround(rgamma(param_dists[aux2], param_dists[1 + aux2]),0);
                  } //end for(l)
                  break;
              } //end switch
            } //end if
            else
              if(cdata[aux_i_1] < 0)
                return;
          } //end for(j)
        } //end for(o)
	if(outfile != NULL)
	{
	  fprintf(of, "i%d.type%d\t", i+1, k+1);
	  for(o = 0 ; o < nro_types ; o++)
	     fprintf(of, "%.0f\t", cdata[o + nro_types * k + nro_types * nro_types * i]);
	  fprintf(of,"\n");
	  fflush(of);
	}
      } //end for(k)
    } //end for(i)

    PutRNGstate();

    if(outfile != NULL)
    {
      fclose(of);
    }
      
  } //end rBGWMindependent

} //end extern
