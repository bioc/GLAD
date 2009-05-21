/*****************************************************************************/
/* Copyright (C) 2009 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2009                            */
/* Contact: glad@curie.fr                                                    */
/*****************************************************************************/

#include "glad-utils.h"

extern "C"
{

  void FilterBkpStep(const int *Chromosome,
		     int *Breakpoints,
		     int *Level,
		     const int *PosOrder,
		     double *NextLogRatio,
		     const double *LogRatio,
		     const int *maxLevel,
		     // ajout de variables pour updateOutliers
		     int *OutliersAws,
		     double *Smoothing,
		     // ajout de variables pour detectOutliers
		     int *OutliersMad,
		     int *OutliersTot,
		     const int *msize,
		     const double *alpha,
		     const int *l,
		     const double *NormalRef,
		     const double *deltaN,
		     int *NormalRange)
  {
  }

  void updateFilterBkp(const int *Chromosome,
		       int *Breakpoints,
		       int *Level,
		       const int *PosOrder,
		       double *NextLogRatio,
		       const double *LogRatio,
		       const int *maxLevel,
		       // ajout de variables pour updateOutliers
		       int *OutliersAws,
		       double *Smoothing,
		       // ajout de variables pour detectOutliers
		       int *OutliersMad,
		       int *OutliersTot,
		       const int *msize,
		       const double *alpha,
		       const int *l,
		       const double *NormalRef,
		       const double *deltaN,
		       int *NormalRange)
  {

    updateLevel(Chromosome,
    		Breakpoints,
    		Level,
    		PosOrder,
    		NextLogRatio,
    		LogRatio,
    		maxLevel,
    		l);


    updateOutliers (OutliersAws,
    		    Level,
    		    Breakpoints,
    		    Smoothing,
    		    l);


    detectOutliers(LogRatio,
    		   Level,
    		   OutliersAws,
    		   OutliersMad,
    		   OutliersTot,
    		   msize,
    		   alpha,
    		   l);

    //recalcul de la smoothing line
    compute_median_smoothing(LogRatio,
			     Level,
			     Smoothing,
			     l);

            
    // on prend comme référence ceux qui sont compris entre certaines valeurs            
    compute_NormalRange(Smoothing,
			NormalRef,
			Level,
			NormalRange,
			deltaN,
			l);



  }
}
