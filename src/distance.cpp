/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998, 2001  Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
 *
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <float.h>
#include <R_ext/Arith.h>
#include <R_ext/Error.h>
#include "mva.h"

extern "C"
{
  double R_euclidean(double *x, int nr, int nc, int i1, int i2)
  {
    double dev, dist;
    int count, j;

    count= 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
      if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	dev = (x[i1] - x[i2]);
	dist += dev * dev;
	count++;
      }
      i1 += nr;
      i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return sqrt(dist);
  }

  double R_maximum(double *x, int nr, int nc, int i1, int i2)
  {
    double dev, dist;
    int count, j;

    count = 0;
    dist = -DBL_MAX;
    for(j = 0 ; j < nc ; j++) {
      if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	dev = fabs(x[i1] - x[i2]);
	if(dev > dist)
	  dist = dev;
	count++;
      }
      i1 += nr;
      i2 += nr;
    }
    if(count == 0) return NA_REAL;
    return dist;
  }

  double R_manhattan(double *x, int nr, int nc, int i1, int i2)
  {
    double dist;
    int count, j;

    count = 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
      if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	dist += fabs(x[i1] - x[i2]);
	count++;
      }
      i1 += nr;
      i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return dist;
  }

  double R_canberra(double *x, int nr, int nc, int i1, int i2)
  {
    double dist, sum, diff;
    int count, j;

    count = 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
      if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	sum = fabs(x[i1] + x[i2]);
	diff = fabs(x[i1] - x[i2]);
	if (sum > DBL_MIN || diff > DBL_MIN) {
	  dist += diff/sum;
	  count++;
	}
      }
      i1 += nr;
      i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return dist;
  }

  double R_dist_binary(double *x, int nr, int nc, int i1, int i2)
  {
    int total, count, dist;
    int j;

    total = 0;
    count = 0;
    dist = 0;

    for(j = 0 ; j < nc ; j++) {
      if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	if(x[i1] || x[i2]){
	  count++;
	  if( ! (x[i1] && x[i2]) ) dist++;
	}
	total++;
      }
      i1 += nr;
      i2 += nr;
    }

    if(total == 0) return NA_REAL;
    if(count == 0) return 0;
    return (double) dist / count;
  }

  /* Pearson / Pearson centered (correlation)
   * Added by A. Lucas
   */

  double R_pearson(double *x, int nr, int nc, int i1, int i2)
  {
    double num,sum1,sum2, dist;
    int count,j;

    count= 0;
    num = 0;
    sum1 = 0;
    sum2 = 0;

    for(j = 0 ; j < nc ; j++) {
      if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	num += (x[i1] * x[i2]);
	sum1 += (x[i1] * x[i1]);
	sum2 += (x[i2] * x[i2]);
	count++;
      }
      i1 += nr;
      i2 += nr;
    }
    if(count == 0) return NA_REAL;
    dist = 1 - ( num / sqrt(sum1 * sum2) );
    return dist;
  }


  double R_correlation(double *x, int nr, int nc, int i1, int i2)
  {
    double num,denum,sumx,sumy,sumxx,sumyy,sumxy;
    int count,j;

    count= 0;
    sumx=0;
    sumy=0;
    sumxx=0;
    sumyy=0;
    sumxy=0;


    for(j = 0 ; j < nc ; j++) {
      if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	sumxy += (x[i1] * x[i2]);
	sumx += x[i1];
	sumy += x[i2];
	sumxx += x[i1] * x[i1];
	sumyy += x[i2] * x[i2];
	count++;
      }
      i1 += nr;
      i2 += nr;
    }
    if(count == 0) return NA_REAL;
    num = sumxy - ( sumx*sumy /count );
    denum = sqrt( (sumxx - (sumx*sumx /count ) )* (sumyy - (sumy*sumy /count ) ) );
    return 1 - (num / denum);
  }

  enum { EUCLIDEAN=1, MAXIMUM, MANHATTAN, CANBERRA, BINARY ,PEARSON, CORRELATION};
  /* == 1,2,..., defined by order in the R function dist */

  void R_distance(double *x, int *nr, int *nc, double *d, int *diag, int *method)
  {
    int dc, i, j, ij;
    double (*distfun)(double*, int, int, int, int) = NULL;

    switch(*method) {
    case EUCLIDEAN:
      distfun = R_euclidean;
      break;
    case MAXIMUM:
      distfun = R_maximum;
      break;
    case MANHATTAN:
      distfun = R_manhattan;
      break;
    case CANBERRA:
      distfun = R_canberra;
      break;
    case BINARY:
      distfun = R_dist_binary;
      break;
    case PEARSON:
      distfun = R_pearson;
      break;
    case CORRELATION:
      distfun = R_correlation;
      break;

    default:
      error("distance(): invalid distance");
    }

    dc = (*diag) ? 0 : 1; /* diag=1:  we do the diagonal */
    ij = 0;
    for(j = 0 ; j <= *nr ; j++)
      for(i = j+dc ; i < *nr ; i++)
	d[ij++] = distfun(x, *nr, *nc, i, j);
  }





}
