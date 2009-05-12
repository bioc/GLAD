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

#ifndef HAARSEG_H_
#define HAARSEG_H_

#include <math.h>
#include <stdlib.h>

extern "C"
{
  int HaarConv(const double * signal, const double * weight, int signalSize,  int stepHalfSize, double * result);
 
  int FindLocalPeaks(const double * signal, int signalSize, int * peakLoc);
 
  int HardThreshold(const double * signal, double threshold, int * peakLoc);

  int UnifyLevels(const int * baseLevel, const int * addonLevel, int windowSize, int signalSize, int * joinedLevel);
 
  int CopyLocVec(const int * source, int * target);

  int AdjustBreaks(const double * signal, int signalSize, const int * peakLoc, int * newPeakLoc);
 
  int StepConv(const double * signal, int signalSize, int pulseSize, double pulseHeight, double * result);

  void rConvAndPeak(const double * signal,
		    const int * signalSize,
		    const int * stepHalfSize,
		    double * convResult,
		    int * peakLoc);

  void rThresAndUnify(const double * addon,
		      const int * signalSize,
		      int * addonPeaks,
		      const int * basePeaks,
		      const double * threshold,
		      const int * windowSize,
		      int * uniPeaks);

  double FDRThres(const double *x, const double q, const double sdev, const int size);

  void SegmentByPeaks(const double *data,
		      const int *peaks,
		      double *segs,
		      const int length_data,
		      const int length_peaks);


  void HaarSegGLAD(const double * signal,
		   const int * signalSize,
		   const int * stepHalfSize,
		   double * convResult,
		   int * peakLoc,
		   const double *breaksFdrQ,
		   const int *haarStartLevel,
		   const int *haarEndLevel,
		   double *segs);
 
}
#endif /*HAARSEG_H_*/
