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
		     int *NormalRange,
		     // paramètres pour findCluster
		     int *ZoneGen,
		     int *method,
		     const double *sigma,
		     const double *d,
		     const double *lambda,
		     const int *nmin,
		     const int *nmax,
		     int *nbclasses,
		     // paramètres pour le calcul du GNL
		     int *ZoneGNL,
		     const double *forceGL1Value,
		     const double *forceGL2Value,
		     const double *NormalRefValue,
		     const double *ampliconValue,
		     const double *deletionValue)


  {
    updateFilterBkp(Chromosome,
		    Breakpoints,
		    Level,
		    PosOrder,
		    NextLogRatio,
		    LogRatio,
		    maxLevel,
		    // ajout de variables pour updateOutliers
		    OutliersAws,
		    Smoothing,
		    // ajout de variables pour detectOutliers
		    OutliersMad,
		    OutliersTot,
		    msize,
		    alpha,
		    l,
		    NormalRef,
		    deltaN,
		    NormalRange);

    findCluster(LogRatio,
		NormalRange,
		OutliersTot,
		ZoneGen,
		method,
		sigma,
		d,
		lambda,
		nmin,
		nmax,
		nbclasses,
		l);


  compute_cluster_LossNormalGain(// variables pour faire la jointure
				      ZoneGen,
				      int *value_dest,
				      l,
				      Smoothing,
				      const double *forceGL1Value,
				      const double *forceGL2Value,
				      const double *NormalRefValue,
				      const double *ampliconValue,
				      const double *deletionValue,
				      //variables pour calcul la médiane par cluster
				      LogRatio,
				      NormalRange);



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
