/* 
 *    Copyright (C) 2008  Erez Ben-Yaacov
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    <http://www.gnu.org/licenses/>
 *
 */

#include <R.h>
#include <vector>
#include <algorithm>
#include <gsl/gsl_cdf.h>

#include "HaarSeg.h"
#include "glad-utils.h"

#define OK 0
#define ERROR_HAARSEG -1
#define NOT_VALID -1
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define abs(X)  ((X) > 0 ? (X) : -(X))
#define MAXINT 2147483647

extern "C"
{
  /*
   * HaarConv : convolve haar wavelet function with a signal, 
   * applying circular padding to the signal.
   * supports weights when weight pointer is not NULL.
   */
  int HaarConv(const double * signal,
	       const double * weight,
	       int signalSize, 
	       int stepHalfSize, 
	       double * result)
  {
    int k;
    int highEnd, lowEnd;
    double stepNorm = 0;
    double lowWeightSum = 0;
    double highWeightSum = 0;
    double lowSquareSum = 0;
    double highSquareSum = 0;
    double lowNonNormed = 0;
    double highNonNormed = 0;
    //    double totalNorm;

    if (stepHalfSize > signalSize) {
      return ERROR_HAARSEG; /* TODO: handle this endcase */
    }
    result[0] = 0;
    if (weight != NULL) {
      /* init weight sums */
      highWeightSum = 0;
      highSquareSum = 0;
      highNonNormed = 0;
      for (k = 0; k < stepHalfSize; k++) {
	highWeightSum += weight[k];
	highSquareSum += weight[k]*weight[k]; 
	highNonNormed += weight[k]*signal[k];
      }
      /* circular padding */
      lowWeightSum = highWeightSum; 
      lowSquareSum = highSquareSum;
      lowNonNormed = -highNonNormed;
    }/*if (weight != NULL) */
    for (k = 1; k < signalSize; k++) {
      highEnd = k + stepHalfSize - 1;
      if (highEnd >= signalSize) {
	highEnd = signalSize - 1 - (highEnd - signalSize);
      }
      lowEnd = k - stepHalfSize - 1;
      if (lowEnd < 0) {
	lowEnd = - lowEnd - 1; 
      }
      if (weight != NULL) {
	lowNonNormed += signal[lowEnd]*weight[lowEnd] - signal[k-1]*weight[k-1];
	highNonNormed += signal[highEnd]*weight[highEnd] - signal[k-1]*weight[k-1]; 
	lowWeightSum += weight[k-1] - weight[lowEnd];
	highWeightSum += weight[highEnd] - weight[k-1];
	lowSquareSum += weight[k-1]*weight[k-1] - weight[lowEnd]*weight[lowEnd];
	highSquareSum += weight[highEnd]*weight[highEnd] - weight[k-1]*weight[k-1];
	result[k] = (lowNonNormed / lowWeightSum + highNonNormed / highWeightSum) * sqrt((double)(stepHalfSize/2));            
	/*            totalNorm = lowSquareSum / (lowWeightSum*lowWeightSum) + highSquareSum / (highWeightSum*highWeightSum); */
	/*            result[k] = (lowNonNormed / lowWeightSum + highNonNormed / highWeightSum) / sqrt(totalNorm); */
      }/*if (weight != NULL) */
      else {
	result[k] = result[k-1] + signal[highEnd] + signal[lowEnd] - 2*signal[k-1];
      }
    }/* for k */
    
    if (weight == NULL) {
      stepNorm = sqrt((double)(2*stepHalfSize));
      for (k = 1; k < signalSize; k++) {
	result[k] /= stepNorm;
      }
    }

    return OK;
  }/* int HaarConv */
 
  /*
   * FindLocalPeaks: find local maxima on positive values,
   * and local minima on negative values.
   * First and last index are never considered extramum.
   */
  int FindLocalPeaks(const double * signal, int signalSize, int * peakLoc)
  {
    int k;
    int maxSuspect, minSuspect;
    int peakLocInd;
	 
    maxSuspect = NOT_VALID;
    minSuspect = NOT_VALID;
    peakLocInd = 0;
    for (k = 1; k < signalSize-1; k++) {
      if (signal[k] > 0) {
	if ((signal[k] > signal[k-1]) && (signal[k] > signal[k+1])) {
	  peakLoc[peakLocInd] = k;
	  peakLocInd++;
	}
	else if ((signal[k] > signal[k-1]) && (signal[k] == signal[k+1])) {
	  maxSuspect = k;
	}
	else if (signal[k] == signal[k-1] && (signal[k] > signal[k+1])) {
	  if (maxSuspect != NOT_VALID) {
	    peakLoc[peakLocInd] = maxSuspect;
	    peakLocInd++;

// 	      for (j = maxSuspect; j <= k; j++) {
// 	      peakLoc[peakLocInd] = j;
// 	      peakLocInd++;
// 	      }
	    maxSuspect = NOT_VALID;					 
	  }
	}
	else if ((signal[k] == signal[k-1]) && (signal[k] < signal[k+1])) {
	  maxSuspect = NOT_VALID;
	}
      }/* if (signal[k] > 0) */
      else if (signal[k] < 0) {
	if ((signal[k] < signal[k-1]) && (signal[k] < signal[k+1])) {
	  peakLoc[peakLocInd] = k;
	  peakLocInd++;
	}
	else if ((signal[k] < signal[k-1]) && (signal[k] == signal[k+1])) {
	  minSuspect = k;
	}
	else if ((signal[k] == signal[k-1]) && (signal[k] < signal[k+1])) {
	  if(minSuspect != NOT_VALID) {
	    peakLoc[peakLocInd] = minSuspect;
	    peakLocInd++;
// 	    /*
// 	      for (j = minSuspect; j <= k; j++) {
// 	      peakLoc[peakLocInd] = j;
// 	      peakLocInd++;
// 	      }/* for j */
	    minSuspect = NOT_VALID;					 
	  }
	}
	else if ((signal[k] == signal[k-1]) && (signal[k] > signal[k+1])) {
	  minSuspect = NOT_VALID;
	}
      }/* else if (signal[k] < 0) */
    }/* for k */
	 
    peakLoc[peakLocInd] = NOT_VALID;
	 
    return OK;
  }/* int FindLocalPeaks */
 
  /*
   * HardThreshold: Apply hard thresholding
   */
  int HardThreshold(const double * signal, 
		    double threshold,
		    int * peakLoc) 
  {
    int k,l;
	 
    k = 0;
    l = 0;
    while (peakLoc[k] != NOT_VALID) {
      if ((signal[peakLoc[k]] >= threshold) || (signal[peakLoc[k]] <= -threshold)) {
	/* peak is over the threshold */
	peakLoc[l] = peakLoc[k];
	l++;
      }
      k++;
    }
    peakLoc[l] = NOT_VALID;
     
    return OK;
  }/* int HardThreshold */
 
  /*
   * UnifyLevels: Unify several decomposition levels
   */
  int UnifyLevels(const int * baseLevel,
		  const int * addonLevel,
		  int windowSize,
		  int signalSize,
		  int * joinedLevel) 
  {
    int baseInd,addonInd,joinedInd;
	 	 
    baseInd = 0;
    addonInd = 0;
    joinedInd = 0;
     
    /* going over all base */
    while (baseLevel[baseInd] != NOT_VALID) {
      while ((addonLevel[addonInd] != NOT_VALID) && 
	     (addonLevel[addonInd] <= (baseLevel[baseInd] + windowSize))) {
	if (addonLevel[addonInd] < (baseLevel[baseInd] - windowSize)) {
	  joinedLevel[joinedInd] = addonLevel[addonInd];
	  joinedInd++;
	}
	addonInd++;
      }/* while ((addonLevel[addonInd] ... */
      joinedLevel[joinedInd] = baseLevel[baseInd];
      joinedInd++;
      baseInd++;
    }/* while (baseLevel[baseInd] */
	 
    /* insert remaining indexes in addon to joined */
    while (addonLevel[addonInd] != NOT_VALID) {
      joinedLevel[joinedInd] = addonLevel[addonInd];
      joinedInd++;
      addonInd++;
    }
    joinedLevel[joinedInd] = NOT_VALID;
	 
    return OK;
  }/* int UnifyLevels */
 
  /*
   * CopyLocVec: copy source index vector to target index vector
   */
  int CopyLocVec(const int * source, int * target) {
    int k;
    k = 0;
    while (source[k] != NOT_VALID) {
      target[k] = source[k];
      k++;
    }
    /*      mexPrintf("CopyLocVec: copied %d elements\n",k); */
    target[k] = NOT_VALID;
     
    return OK;
  }/* int CopyLocVec */
 
  /*
   *  AdjustBreaks: improving localization of breaks by using a suboptimal,
   *  linear complexity procedure. We try to move each break 1 sample 
   *  left/right, choosing the offset which leads to minimum data error. 	
   */
  int AdjustBreaks(const double * signal,
		   int signalSize,
		   const int * peakLoc,
		   int * newPeakLoc) {
    int k,m,p;
    int n1,n2;
    int bestOffset;
    double s1, s2, ss1, ss2;
    double score, bestScore;
    
    k = 0;
    while (peakLoc[k] != NOT_VALID) {
      newPeakLoc[k] = peakLoc[k];
      k++;
    }
    newPeakLoc[k] = NOT_VALID;
    
    k = 0;
    n1 = 0;
    n2 = 0;
    while (newPeakLoc[k] != NOT_VALID) {
      /* calculating width of segments around the breakpoint */
      if (k == 0) {
	n1 = newPeakLoc[k];
      }
      else {
	n1 = newPeakLoc[k] - newPeakLoc[k-1];
      }
      if (newPeakLoc[k+1] == NOT_VALID) {
	n2 = signalSize - newPeakLoc[k];
      }
      else {
	n2 = newPeakLoc[k+1] - newPeakLoc[k];
      }
        
      /* finding the best offset for current breakpoint, trying only 1 sample offset */
      bestScore = MAXINT;
      bestOffset = 0;
      for (p = -1; p <= 1; p++) {
	/* pointless to try and remove single sample segments */
	if ((n1 == 1) && (p == -1)) {
	  continue;
	}
	if ((n2 == 1) && (p == 1)) {
	  continue;
	}
                
	s1 = 0;
	for (m = (newPeakLoc[k] - n1); m <= (newPeakLoc[k] + p - 1); m++) {
	  s1 += signal[m];
	}
	s1 = s1 / (n1 + p);
	s2 = 0;
	for (m = (newPeakLoc[k] + p); m <= (newPeakLoc[k] + n2 - 1); m++) {
	  s2 += signal[m];
	}
	s2 = s2 / (n2 - p);
            
	ss1 = 0;
	for (m = (newPeakLoc[k] - n1); m <= (newPeakLoc[k] + p - 1); m++) {
	  ss1 += (signal[m] - s1)*(signal[m] - s1);
	}
	ss2 = 0;
	for (m = (newPeakLoc[k] + p); m <= (newPeakLoc[k] + n2 - 1); m++) {
	  ss2 += (signal[m] - s2)*(signal[m] - s2);
	}
	score = ss1 + ss2;
	if (score < bestScore) {
	  bestScore = score;
	  bestOffset = p;
	}
      }/* for p */
      newPeakLoc[k] += bestOffset; 
      k++;
    }/* while newPeakLoc */
    
    return OK;
  }/* int AdjustBreaks */
 
  /*
   * StepConv : convolve a pulse function with a signal, 
   * applying circular padding to the signal.
   */
  int PulseConv(const double * signal,
		int signalSize, 
		int pulseSize,
		double pulseHeight, 
		double * result)
  {
    int k, n, tail, head;

    if (pulseSize > signalSize) {
      return ERROR_HAARSEG; /* TODO: handle this endcase */
    }
    /* circular padding init */
    result[0] = 0;
    for (k = 0; k < ((pulseSize + 1)/2); k++) {
      result[0] += signal[k];
    }
    for (k = 0; k < (pulseSize/2); k++) {
      result[0] += signal[k];
    }  
    result[0] *= pulseHeight;  
    n = 1;
    for (k = (pulseSize/2); k < signalSize + (pulseSize/2) - 1; k++) {
      tail = k - pulseSize;
      if (tail < 0) {
	tail = -tail - 1; 
      }
      head = k;
      if (head >= signalSize) {
	head = signalSize - 1 - (head - signalSize);
      }
      result[n] = result[n-1] + ((signal[head] - signal[tail]) * pulseHeight);
      n++;
    }/* for k */
    
    return OK;
  }/* int PulseConv */                 

  /*************************/
  /***  HaarSegGLAD   ******/
  /*************************/

  void HaarSegGLAD(const double * signal,
		   const int * signalSize,
		   const int * stepHalfSize,
		   double * convResult,
		   int * peakLoc,
		   const double *breaksFdrQ,
		   const int *haarStartLevel,
		   const int *haarEndLevel,
		   double *segs,
		   const double *weights)

  {
    double *convResult_tmp1;
    double *fdr_tmp1;
    double peakSigmaEst;
    double T;

    int unifyWin;
    int level;
    int stepHalfSize_tmp1;
    int i;
    int indice;

    int *uniPeakLoc;
    int *tmpPeakLoc;
    int *breakpoints;
    int *peakLoc_tmp1;

    const int size = *signalSize;


    convResult_tmp1 = (double *)calloc(size, sizeof(double));
    peakLoc_tmp1 = (int *)calloc(size, sizeof(double));
    uniPeakLoc = (int *)calloc(size, sizeof(double));
    tmpPeakLoc = (int *)calloc(size, sizeof(double));

    uniPeakLoc[0] = -1;

    rConvAndPeak(signal,
		 &size,
		 stepHalfSize,
		 convResult,
		 peakLoc);

    peakSigmaEst = median_fabs_double(convResult, size) / 0.6745;

    for (level = *haarStartLevel; level <= *haarEndLevel; level++)
      {
	stepHalfSize_tmp1 = (int) pow(2, (double)level);

	if(weights == NULL)
	  {
	    printf("Il n'y a pas de poids\n");
	    rConvAndPeak(signal,
			 &size,
			 &stepHalfSize_tmp1,
			 convResult_tmp1,
			 peakLoc_tmp1);

	  }
	else
	  {
	    rWConvAndPeak(signal,
			  weights,
			  &size,
			  &stepHalfSize_tmp1,
			  convResult_tmp1,
			  peakLoc_tmp1);

	    printf("WWWWWWWWWWWWWWWWW\n");
	  }



	for (i = 0; i < size; i++)
	  {
	    if(peakLoc_tmp1[i] == -1)
	      {
		break;
	      }
	  }

	indice = --i;

	T = 0;
	if (indice >= 0)
	  {
	    fdr_tmp1 = (double *)malloc((indice + 1) * sizeof(double));
	    for(i = 0; i < (indice + 1); i++)
	      {
		fdr_tmp1[i] = convResult_tmp1[peakLoc_tmp1[i]];
	      }

	    T = FDRThres(fdr_tmp1, *breaksFdrQ, peakSigmaEst, indice + 1);
	    free(fdr_tmp1);
	  }


	unifyWin = (int) pow(2, (double)(level - 1));

	memcpy(tmpPeakLoc, uniPeakLoc, size * sizeof(int)); 
	
	for (i = 0; i < size; i++)
	  {
	    uniPeakLoc[i] = 0;
	  }

	rThresAndUnify(convResult_tmp1, 
		       &size, 
		       peakLoc_tmp1,
		       tmpPeakLoc,
		       &T,
		       &unifyWin,
		       uniPeakLoc);

      } // fin de la boucle sur les levels

    for (i = 0; i < size; i++)
      {
	if(uniPeakLoc[i] == -1)
	  {
	    break;
	  }
      }

    indice = --i;

    breakpoints = (int *)calloc((indice + 1), sizeof(int));
    if (indice >= 0)
      {
 	for(i = 0; i < (indice + 1); i++)
 	  {
	    breakpoints[i] = uniPeakLoc[i]; // on ne fait pas "+1" pour avoir des indices qui commencent à 0 dans SegmentByPeaks
 	  }
      }


    SegmentByPeaks(signal,
		   breakpoints,
		   segs,
		   size,
		   indice + 1);


    free(breakpoints);
    free(convResult_tmp1);
    free(peakLoc_tmp1);
    free(uniPeakLoc);
    free(tmpPeakLoc);

  }//rConvAndPeak



  bool plusgrand (double i,double j) 
  { 
    return (i > j); 
  }

  double FDRThres(const double *x, const double q, const double sdev, const int size)
  {
    int i, k = -1;
    double m;
    double p;
    vector<double> sortedX;


    if (size < 2)
      {
	return 0;
      }
    else
      {
	for(i = 0; i < size; i++)
	  {
	    sortedX.push_back(fabs(x[i]));
	  }

	// On trie le vecteur vec
	sort(sortedX.begin(), sortedX.end(), plusgrand);
      }

    for(i = 0; i < size; i++)
      {
	m = (double)(i + 1) / size;
	p = 2 - 2 * gsl_cdf_gaussian_P (sortedX[i], sdev);

	if (p <= (m * q))
	  {
	    k = i;
	  }
      }

    if(k == -1)
      {
	return sortedX[0] + 0.0000000000000001;
      }
    else
      {
	return sortedX[k];
      }
  }

  void SegmentByPeaks(const double *data,
		      const int *peaks,
		      double *segs,
		      const int length_data,
		      const int length_peaks)
  {

    // pour l'instant cette fonction ne prend pas en compte les poids
    int i;
    int k;
    int *st;
    int *ed;

    double sum;

    st = (int *)malloc((length_peaks + 1) * sizeof(int));
    st[0] = 0;
    memcpy(&st[1], peaks, length_peaks * sizeof(int));

    ed = (int *)malloc((length_peaks + 1) * sizeof(int));
    ed[length_peaks] = length_data - 1;

    for(i = 0; i < length_peaks; i++)
      {
	ed[i] = peaks[i] - 1;
      }


    for (k = 0; k <= length_peaks; k++ )
      {
	sum = 0;
	for (i = st[k]; i <= ed[k]; i++)
	  {
	    sum += data[i];
	  }

	sum /= (double) (ed[k] - st[k] + 1);


	for (i = st[k]; i <= ed[k]; i++)
	  {
	    segs[i] = sum;
	  }
      }


    free(st);
    free(ed);

  }


  /*****************************************/
  /*** contenu du fichier t_haarseg   ******/
  /*****************************************/

  void rConvAndPeak(const double * signal,
		    const int * signalSize,
		    const int * stepHalfSize,
		    double * convResult,
		    int * peakLoc) 
  {
    HaarConv(signal, NULL, *signalSize, *stepHalfSize, convResult);
    FindLocalPeaks(convResult, *signalSize, peakLoc);                     
  }//rConvAndPeak

  void rWConvAndPeak(const double * signal,
		     const double * weight,
		     const int * signalSize,
		     const int * stepHalfSize,
		     double * convResult,
		     int * peakLoc) {
    HaarConv(signal, weight, *signalSize, *stepHalfSize, convResult);
    FindLocalPeaks(convResult, *signalSize, peakLoc);                     
  }//rWConvAndPeak


  void rThresAndUnify(const double * addon,
		      const int * signalSize,
		      int * addonPeaks,
		      const int * basePeaks,
		      const double * threshold,
		      const int * windowSize,
		      int * uniPeaks) {
    HardThreshold(addon, *threshold, addonPeaks);
    UnifyLevels(basePeaks, addonPeaks, *windowSize, *signalSize, uniPeaks);
  }//rThresAndUnify

  void rAdjustBreaks(const double * signal,
		     const int * signalSize,
		     const int * peakLoc,
		     int * newPeakLoc) {
    AdjustBreaks(signal, *signalSize, peakLoc, newPeakLoc);
  }//rAdjustBreaks

  void rPulseConv(const double * signal,
		  const int * signalSize, 
		  const int * pulseSize,
		  const double * pulseHeight, 
		  double * result) {
    PulseConv(signal, *signalSize, *pulseSize, *pulseHeight, result);
  }//rPulseConv

}
