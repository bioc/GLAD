/*
 *  R : A Computer Language for Statistical Data Analysis

 *  Copyright (C) 1999-2008   The R Development Core Team
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
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/.
 */

#include <stdio.h>
#include <stdlib.h>

// #include <R_ext/Boolean.h>
// #include <Rinternals.h>
// #include "mva.h"

// SEXP R_cutree(SEXP merge, SEXP which)
// {
//   /* Return grouping vector from cutting a (binary) (cluster) tree
//    * into which[j] groups.
//    * merge = (n-1) x 2  matrix, described in help(hclust) */
//   SEXP ans;
//   int n, k, l, nclust, m1, m2, j, mm = 0;
//   Rboolean found_j, *sing;
//   int *m_nr, *z;

//   merge = coerceVector(merge, INTSXP);
//   which = coerceVector(which, INTSXP);

//   n = nrows(merge)+1;
//   /* using 1-based indices ==> "--" */
//   sing = (Rboolean *) R_alloc(n, sizeof(Rboolean)); sing--;
//   m_nr = (int *) R_alloc(n, sizeof(int)); m_nr--;
//   z	 = (int *) R_alloc(n, sizeof(int)); z--;
//   PROTECT(ans = allocMatrix(INTSXP, n, LENGTH(which)));

//   for(k = 1; k <= n; k++) {
//     sing[k] = TRUE;/* is k-th obs. still alone in cluster ? */
//     m_nr[k] = 0;/* containing last merge-step number of k-th obs. */
//   }

//   for(k = 1; k <= n-1; k++) {
//     /* k-th merge, from n-k+1 to n-k atoms: (m1,m2) = merge[ k , ] */
//     m1 = INTEGER(merge)[k-1];
//     m2 = INTEGER(merge)[n-1+k-1];

//     if(m1 < 0 && m2 < 0) {/* merging atoms [-m1] and [-m2] */
//       m_nr[-m1] = m_nr[-m2] = k;
//       sing[-m1] = sing[-m2] = FALSE;
//     }
//     else if(m1 < 0 || m2 < 0) {/* the other >= 0 */
//       if(m1 < 0) { j = -m1; m1 = m2; } else j = -m2;
//       /* merging atom j & cluster m1 */
//       for(l = 1; l <= n; l++)
// 	if (m_nr[l] == m1)
// 	  m_nr[l] = k;
//       m_nr[j] = k;
//       sing[j] = FALSE;
//     }
//     else { /* both m1, m2 >= 0 */
//       for(l=1; l <= n; l++) {
// 	if( m_nr[l] == m1 || m_nr[l] == m2 )
// 	  m_nr[l] = k;
//       }
//     }

//     /* does this k-th merge belong to a desired group size which[j] ?
//      * if yes, find j (maybe multiple ones): */
//     found_j = FALSE;
//     for(j = 0; j < LENGTH(which); j++) {
//       if(INTEGER(which)[j] == n - k) {
// 	if(!found_j) { /* first match (and usually only one) */
// 	  found_j = TRUE;
// 	  for(l = 1; l <= n; l++)
// 	    z[l] = 0;
// 	  nclust = 0;
// 	  mm = j*n; /*may want to copy this column of ans[] */
// 	  for(l = 1, m1 = mm; l <= n; l++, m1++) {
// 	    if(sing[l])
// 	      INTEGER(ans)[m1] = ++nclust;
// 	    else {
// 	      if (z[m_nr[l]] == 0)
// 		z[m_nr[l]] = ++nclust;
// 	      INTEGER(ans)[m1] = z[m_nr[l]];
// 	    }
// 	  }
// 	}
// 	else { /* found_j: another which[j] == n-k : copy column */
// 	  for(l = 1, m1 = j*n, m2 = mm; l <= n; l++, m1++, m2++)
// 	    INTEGER(ans)[m1] = INTEGER(ans)[m2];
// 	}
//       } /* if ( match ) */
//     } /* for(j .. which[j] ) */
//   } /* for(k ..) {merge} */

//     /* Dealing with trivial case which[] = n : */
//   for(j = 0; j < LENGTH(which); j++)
//     if(INTEGER(which)[j] == n)
//       for(l = 1, m1 = j*n; l <= n; l++, m1++)
// 	INTEGER(ans)[m1] = l;

//   UNPROTECT(1);
//   return(ans);
// }


extern "C" 
{
  void R_cutree(int *merge, int *which, int *ans, int *nbelt)
  {

    // n est le nombre déléments
    // which est le nombre de classes que l'on veut

    /* Return grouping vector from cutting a (binary) (cluster) tree
     * into which[j] groups.
     * merge = (n-1) x 2  matrix, described in help(hclust) */
    //    SEXP ans;

    const int n = *nbelt;
    int  k, l, nclust, m1, m2, j, mm = 0;
    int found_j, *sing;
    int *m_nr, *z;

    //    merge = coerceVector(merge, INTSXP);
    //    which = coerceVector(which, INTSXP);

    //    n = nrows(merge)+1;
    /* using 1-based indices ==> "--" */
    //     sing = (Rboolean *) R_alloc(n, sizeof(Rboolean)); sing--;
    //     m_nr = (int *) R_alloc(n, sizeof(int)); m_nr--;
    //     z	 = (int *) R_alloc(n, sizeof(int)); z--;
    sing = (int *)calloc(n, sizeof(int));
    sing--;
    m_nr = (int *)calloc(n, sizeof(int));
    m_nr--;
    z = (int *)calloc(n, sizeof(int));
    z--;
    //    PROTECT(ans = allocMatrix(INTSXP, n, LENGTH(which)));

    for(k = 1; k <= n; k++) 
      {
	sing[k] = 1;/* is k-th obs. still alone in cluster ? */
	m_nr[k] = 0;/* containing last merge-step number of k-th obs. */
      }

    for(k = 1; k <= n-1; k++) 
      {
	/* k-th merge, from n-k+1 to n-k atoms: (m1,m2) = merge[ k , ] */
	m1 = merge[k-1];
	m2 = merge[n-1+k-1];

	if(m1 < 0 && m2 < 0) {/* merging atoms [-m1] and [-m2] */
	  m_nr[-m1] = m_nr[-m2] = k;
	  sing[-m1] = sing[-m2] = 0;
	}
	else if(m1 < 0 || m2 < 0) {/* the other >= 0 */
	  if(m1 < 0) { j = -m1; m1 = m2; } else j = -m2;
	  /* merging atom j & cluster m1 */
	  for(l = 1; l <= n; l++)
	    if (m_nr[l] == m1)
	      m_nr[l] = k;
	  m_nr[j] = k;
	  sing[j] = 0;
	}
	else { /* both m1, m2 >= 0 */
	  for(l=1; l <= n; l++) {
	    if( m_nr[l] == m1 || m_nr[l] == m2 )
	      m_nr[l] = k;
	  }
	}

      /* does this k-th merge belong to a desired group size which[j] ?
       * if yes, find j (maybe multiple ones): */
      found_j = 0;
      for(j = 0; j < 1; j++) {
	if(which[j] == n - k) {
	  if(!found_j) { /* first match (and usually only one) */
	    found_j = 1;
	    for(l = 1; l <= n; l++)
	      z[l] = 0;
	    nclust = 0;
	    mm = j*n; /*may want to copy this column of ans[] */
	    for(l = 1, m1 = mm; l <= n; l++, m1++) {
	      if(sing[l])
		ans[m1] = ++nclust;
	      else {
		if (z[m_nr[l]] == 0)
		  z[m_nr[l]] = ++nclust;
		ans[m1] = z[m_nr[l]];
	      }
	    }
	  }
	  else { /* found_j: another which[j] == n-k : copy column */
	    for(l = 1, m1 = j*n, m2 = mm; l <= n; l++, m1++, m2++)
	      ans[m1] = ans[m2];
	  }
	} /* if ( match ) */
      } /* for(j .. which[j] ) */
    } /* for(k ..) {merge} */

    /* Dealing with trivial case which[] = n : */
    for(j = 0; j < 1; j++)
      if(which[j] == n)
	for(l = 1, m1 = j*n; l <= n; l++, m1++)
	  ans[m1] = l;

    sing++;
    m_nr++;
    z++;
    free(sing);
    free(m_nr);
    free(z);

  }
}
