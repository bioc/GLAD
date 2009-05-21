/*****************************************************************************/
/* Copyright (C) 2009 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2009                            */
/* Contact: glad@curie.fr                                                    */
/*****************************************************************************/

#include "loopRemove.h"
#include "findCluster.h"
#include "glad-utils.h"

extern "C"
{
  void daglad_OptmisationBreakpoints_findCluster(double *Smoothing,
						 int *NormalRange,
						 const double *NormalRef,
						 const double *deltaN,
						 // variable pour loop_chromosome_removeLevel
						 const double *LogRatio,
						 double *NextLogRatio,
						 const int *PosOrder,
						 int *Level,
						 int *OutliersAws,
						 int *OutliersMad,
						 int *OutliersTot,
						 int *Breakpoints,
						 const int *msize,
						 const double *alpha,
						 const double *lambda,
						 const double *d,
						 const double *sigma,
						 const int *NbChr,   // Nombre de chromosome à analyser
						 const int *sizeChr, // taille de chaque chromosome
						 const int *startChr, // position pour le début des valeurs de chaque chromosome
						 const int *BkpDetected,
						 // paramètres pour findCluster
						 const int *method, // méthode de clustering
						 const double *findClusterSigma,
						 const double *lambdaclusterGen,
						 const int *nmin,
						 const int *nmax,
						 int *ZoneGen,
						 int *nbclasses,
						 // paramètres pour le calcul du GNL
						 int *ZoneGNL,
						 const double *forceGL1Value,
						 const double *forceGL2Value,
						 const double *NormalRefValue,
						 const double *ampliconValue,
						 const double *deletionValue,
						 const int *l) // nombre total de sondes
  {

    OptmisationBreakpointsStep(Smoothing,
			       NormalRange,
			       NormalRef,
			       deltaN,
			       // variable pour loop_chromosome_removeLevel
			       LogRatio,
			       NextLogRatio,
			       PosOrder,
			       Level,
			       OutliersAws,
			       OutliersMad,
			       OutliersTot,
			       Breakpoints,
			       msize,
			       alpha,
			       lambda,
			       d,
			       sigma,
			       NbChr,   // Nombre de chromosome à analyser
			       sizeChr, // taille de chaque chromosome
			       startChr, // position pour le début des valeurs de chaque chromosome
			       BkpDetected,
			       l); // nombre total de sondes


    findCluster(LogRatio,
		NormalRange,
		OutliersTot,
		ZoneGen,
		method,
		// paramètres pour clusterglad
		findClusterSigma,
		d,
		lambdaclusterGen,
		nmin,
		nmax,
		nbclasses,
		l);


  void compute_cluster_LossNormalGain(// variables pour faire la jointure
				      ZoneGen,
				      ZoneGNL,
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
}
    
