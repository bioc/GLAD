/*****************************************************************************/
/* Copyright (C) 2009 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2009                            */
/* Contact: glad@curie.fr                                                    */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#ifdef IS_MAC_OS
#include <limits.h>
#else
#include <values.h>
#endif

#ifndef MAXINT
#define MAXINT INT_MAX
#endif

#include "MoveBkp.h"
#include "findCluster.h"
#include "glad-utils.h"

extern "C"
{

  void MoveBkp_updateOutliers (int *OutliersAws,
			       int *OutliersTot,
			       int *Level,
			       int *Region,
			       int *Breakpoints,
			       double *Smoothing,
			       int *ZoneGNL,
			       const int *l)
  {
  
    int pos;
    int pos_moins_un;
    int last_bkp_pos = -1;
    const int nb = *l - 1;

    for (pos = 1; pos < nb; pos++)
      {
	pos_moins_un = pos-1;      
	if (Level[pos_moins_un] == Level[pos+1] && 
	    Level[pos_moins_un] != Level[pos])
	  {
	    Level[pos] = Level[pos_moins_un];
	    Region[pos] = Region[pos_moins_un];
	    Breakpoints[pos_moins_un] = 0;
	    Breakpoints[pos] = 0;
	    OutliersAws[pos] = 1;
	    OutliersTot[pos] = 1;
	    ZoneGNL[pos] = ZoneGNL[pos_moins_un];
	    Smoothing[pos] = Smoothing[pos_moins_un];
	  }

	if(Breakpoints[pos] == 1)
	  {
	    if(pos_moins_un == last_bkp_pos)
	      {
		if(pos > 1)
		  {
		    if (Smoothing[pos] > Smoothing[pos_moins_un])
		      {            
			OutliersAws[pos] = 1;
			OutliersTot[pos] = 1;
			Level[pos] = Level[pos_moins_un];
			Region[pos] = Region[pos_moins_un];
		      }
		    else
		      {
			OutliersAws[pos] = -1;
			OutliersTot[pos] = -1;
			Level[pos] = Level[pos_moins_un];
			Region[pos] = Region[pos_moins_un];
		      }
		    Breakpoints[pos_moins_un] = 0;
		  }
	      }
	    last_bkp_pos=pos;
	  }
      }
  }




  void MoveBkp_StepAll(const int *subBkpInfo_MoveBkp,
		       const int *subBkpInfo_PosOrder,
		       const double *LogRatio,
		       double *NextLogRatio,
		       const int *Chromosome,
		       const int *PosOrder,
		       int *Breakpoints,
		       int *OutliersTot,
		       int *OutliersAws,
		       int *OutliersMad,
		       int *Level,
		       int *Region,
		       double *Smoothing,
		       int *ZoneGNL,
		       int *NormalRange,
		       // seuils
		       const double *NormalRef,
		       const double *deltaN,
		       const double *forceGL1Value,
		       const double *forceGL2Value,
		       const double *ampliconValue,
		       const double *deletionValue,
		       // paramètres pour findCluster
		       int *method,
		       const double *Sigma,
		       const double *d,
		       const double *lambda,
		       const int *min,
		       const int *max,
		       int *nbclasses,
		       const int *lensubBkp,
		       const int *l)
  {

    int *zone;

    zone = (int *)malloc(*l * sizeof(int));

    MoveBkp_Step1(subBkpInfo_MoveBkp,
		  subBkpInfo_PosOrder,
		  LogRatio,
		  NextLogRatio,
		  Chromosome,
		  PosOrder,
		  Breakpoints,
		  OutliersTot,
		  OutliersAws,
		  OutliersMad,
		  Level,
		  Region,
		  Smoothing,
		  ZoneGNL,
		  NormalRange,
		  NormalRef,
		  deltaN,
		  lensubBkp,
		  l);


    // on fait le clustering sur NormalRange
    findCluster(LogRatio,
		NormalRange,              
		OutliersTot,
		zone,
		method,
		// paramètres pour clusterglad
		Sigma,
		d,
		lambda,
		min,
		max,
		nbclasses,
		l);


    MoveBkp_Step2(OutliersAws,
		  OutliersTot,
		  Level,
		  Region,
		  Breakpoints,
		  // variables pour faire la jointure
		  ZoneGNL,
		  zone,
		  l,
		  Smoothing,
		  forceGL1Value,
		  forceGL2Value,
		  NormalRef,
		  ampliconValue,
		  deletionValue,
		  deltaN,
		  //variables pour calcul la médiane par cluster
		  LogRatio,
		  NormalRange);


    // on fait le clustering sur NormalRange
    findCluster(LogRatio,
		NormalRange,              
		OutliersTot,
		zone,
		method,
		// paramètres pour clusterglad
		Sigma,
		d,
		lambda,
		min,
		max,
		nbclasses,
		l);


    compute_cluster_LossNormalGain(// variables pour la jointure
				   zone,
				   ZoneGNL,
				   l,
				   Smoothing,
				   forceGL1Value,
				   forceGL2Value,
				   NormalRef,
				   ampliconValue,
				   deletionValue,  
				   // variables pour le calcul de la médiane par cluster
				   LogRatio,
				   NormalRange);



    free(zone);

  }

  void MoveBkp_Step1(const int *subBkpInfo_MoveBkp,
		     const int *subBkpInfo_PosOrder,
		     const double *LogRatio,
		     double *NextLogRatio,
		     const int *Chromosome,
		     const int *PosOrder,
		     int *Breakpoints,
		     int *OutliersTot,
		     int *OutliersAws,
		     int *OutliersMad,
		     int *Level,
		     int *Region,
		     double *Smoothing,
		     int *GNL,
		     int *NormalRange,
		     const double *NormalRef,
		     const double *deltaN,
		     const int *lensubBkp,
		     const int *l)
  {

    int i;
    int max_Level = -MAXINT;

    MoveBkp_Delete_Bkp(subBkpInfo_MoveBkp,
		       subBkpInfo_PosOrder,
		       Breakpoints,
		       OutliersTot,
		       OutliersAws,
		       OutliersMad,
		       Level,
		       Region,
		       Smoothing,
		       GNL,
		       lensubBkp);


    // recalcul de la smoothing line
    compute_median_smoothing(LogRatio,
			     Level,
			     Smoothing,
			     l);


    for (i = 0; i < *l; i++)
      {
	if(Level[i] > max_Level)
	  {
	    max_Level = Level[i];
	  }
      }


    updateLevel(Chromosome,
		Breakpoints,
		Level,
		PosOrder,
		NextLogRatio,
		LogRatio,
		&max_Level,
		l);


    compute_NormalRange(Smoothing,
			NormalRef,
			Level,
			NormalRange,
			deltaN,
			l);
                           


  }


  void MoveBkp_Step2(int *OutliersAws,
		     int *OutliersTot,
		     int *Level,
		     int *Region,
		     int *Breakpoints,
		     // variables pour faire la jointure
		     int *ZoneGNL,
		     int *value_dest,
		     const int *length_dest,
		     double *Smoothing,
		     const double *forceGL1Value,
		     const double *forceGL2Value,
		     const double *NormalRefValue,
		     const double *ampliconValue,
		     const double *deletionValue,
		     const double *deltaN,
		     //variables pour calcul la médiane par cluster
		     const double *LogRatio,
		     int *NormalRange)
  {

    const int *l = length_dest;


    compute_cluster_LossNormalGain(// variables pour faire la jointure
				   ZoneGNL,
				   value_dest,
				   length_dest,
				   Smoothing,
				   forceGL1Value,
				   forceGL2Value,
				   NormalRefValue,
				   ampliconValue,
				   deletionValue,
				   //variables pour calcul la médiane par cluster
				   LogRatio,
				   NormalRange);

    // Mise à jour des outliers
    MoveBkp_updateOutliers(OutliersAws,
			   OutliersTot,
			   Level,
			   Region,
			   Breakpoints,
			   Smoothing,
			   ZoneGNL,
			   l);

    // recalcul de la smoothing line
    compute_median_smoothing(LogRatio,
			     Level,
			     Smoothing,
			     l);



    compute_NormalRange(Smoothing,
			NormalRefValue,
			Level,
			NormalRange,
			deltaN,
			l);
  }


  void MoveBkp_Delete_Bkp(const int *subBkpInfo_MoveBkp,
			  const int *subBkpInfo_PosOrder,
			  int *Breakpoints,
			  int *OutliersTot,
			  int *OutliersAws,
			  int *OutliersMad,
			  int *Level,
			  int *Region,
			  double *Smoothing,
			  int *GNL,
			  const int *l)
  {

    int i;
    int indexPos;
    int indexPosNext;
    int indexPosBefore;
    const int nb = *l;

    for (i = 0; i < nb; i++)
      {
	if (subBkpInfo_MoveBkp[i] == 1)
	  {
	    // déplacement à droite
	    indexPos = subBkpInfo_PosOrder[i] - 1;
	    indexPosNext = indexPos + 1;
	    Breakpoints[indexPos] = 0;
	    Breakpoints[indexPosNext] = 1;
	    OutliersTot[indexPosNext] = 0;
	    OutliersAws[indexPosNext] = 0;
	    OutliersMad[indexPosNext] = 0;
	    OutliersTot[indexPos] = 0;
	    OutliersAws[indexPos] = 0;
	    OutliersMad[indexPos] = 0;
	    Level[indexPosNext] = Level[indexPos];
	    Region[indexPosNext] = Region[indexPos];
	    Smoothing[indexPosNext] = Smoothing[indexPos];
	    GNL[indexPosNext] = GNL[indexPos];
	  }
	else
	  {
	    if (subBkpInfo_MoveBkp[i] == -1)
	      {
		// déplacement à gauche
		indexPos = subBkpInfo_PosOrder[i] - 1;
		indexPosBefore = indexPos - 1;
		indexPosNext = indexPos + 1;
		Breakpoints[indexPos] = 0;
		Breakpoints[indexPosBefore] = 1;
		OutliersTot[indexPosBefore] = 0;
		OutliersAws[indexPosBefore] = 0;
		OutliersMad[indexPosBefore] = 0;
		OutliersTot[indexPos] = 0;
		OutliersAws[indexPos] = 0;
		OutliersMad[indexPos] = 0;
		Level[indexPos] = Level[indexPosNext];
		Region[indexPosNext] = Region[indexPos];
		Smoothing[indexPos] = Smoothing[indexPosNext];
		GNL[indexPos] = GNL[indexPosNext];
	      }
	  }
      }
  }

  void MoveBkp_updateGNL(int *ZoneGNL,
			 const double *Smoothing,
			 const int *OutliersTot,
			 const int *l)
  {
    double *minG;
    double *maxL;
    double *minAmp;
    double *maxDel;
    int i;
    const int nb=*l;

    minG=(double *)malloc(1 * sizeof(double));
    maxL=(double *)malloc(1 * sizeof(double));
    minAmp=(double *)malloc(1 * sizeof(double));
    maxDel=(double *)malloc(1 * sizeof(double));

    rangeGainLoss(Smoothing, ZoneGNL, OutliersTot, minG, maxL, minAmp, maxDel, l);


    for (i = 0; i < nb; i++)
      {
	ZoneGNL[i] = 0;

	if (Smoothing[i] >= *minG || Smoothing[i] >= *minAmp)
	  {
	    ZoneGNL[i] = 1;
	 
	    if(Smoothing[i] >= *minAmp)
	      {
		ZoneGNL[i] = 2;
	      }

	  }

	else
	  {
	    if (Smoothing[i] <= *maxL || Smoothing[i] <= *maxDel)
	      {
		ZoneGNL[i]=-1;

		if(Smoothing[i] <= *maxDel)
		  {
		    ZoneGNL[i] = -10;
		  }
	      }
	  }
      }

    free(minG);
    free(maxL);
    free(minAmp);
    free(maxDel);

  }

  void rangeGainLoss(const double *Smoothing,
		     const int *ZoneGNL,
		     const int *OutliersTot,
		     double *minG,
		     double *maxL,
		     double *minAmp,
		     double *maxDel,
		     const int *l)
  {
    int i;
    const int nb = *l;

    *minG = MAXDOUBLE;
    *minAmp = MAXDOUBLE;
    *maxL = -MAXDOUBLE;
    *maxDel = -MAXDOUBLE;

    for (i = 0; i < nb; i++)
      {
	if (OutliersTot[i] == 0)
	  {

	    switch(ZoneGNL[i])
	      {
	      case 0:
		break;

	      case 1:
		if(Smoothing[i] < *minG)
		  {
		    *minG = Smoothing[i];
		  }
		break;

	      case -1:
		if (Smoothing[i] > *maxL)
		  { 
		    *maxL = Smoothing[i]; 
		  }
		break;

	      case 2:
		if(Smoothing[i] < *minAmp)
		  {
		    *minAmp = Smoothing[i];
		  }
		break;

	      case -10:
		if(Smoothing[i] > *maxDel)
		  {
		    *maxDel = Smoothing[i];
		  }
		break;

	      }
	  }

      }
  }


  /*******************************************/
  /* fonctions utilisées dans testMoveBkp.R  */
  /*******************************************/

  void loopTestBkpToMove(const double *LogRatio,
			 const double *NextLogRatio,
			 const double *Smoothing,
			 const double *SmoothingNext,
			 const int *PosOrder,
			 const int *MaxPosOrder,
			 const int *MinPosOrder,
			 int *MoveBkp,
			 const int *NbBkp)
  {
    int i;
    int l=*NbBkp;

    for (i=0;i<l;i++)
      {
	MoveBkp[i]=testSingle(LogRatio[i],NextLogRatio[i],Smoothing[i],SmoothingNext[i]);

	if(MoveBkp[i]==1)
	  {
	    if((PosOrder[i]+1)==MaxPosOrder[i])
	      {
		MoveBkp[i]=0;
	      }
	  }
	else
	  {
	    if(MoveBkp[i]==-1)
	      {
		if((PosOrder[i]-1)==MinPosOrder[i])
		  {
		    MoveBkp[i]=0;
		  }
	      }
	  }
      }
  }



  int  testSingle(const double LogRatio,
		  const double NextLogRatio,
		  const double Smoothing,
		  const double SmoothingNext)
  {
    int moveBkp=0;
    // le créneau est plus bas à droite qu'à gauche
    if (Smoothing > SmoothingNext)
      {
	if ((SmoothingNext <= LogRatio) && (LogRatio <= Smoothing))
	  {
	    if((LogRatio-SmoothingNext) < (Smoothing-LogRatio))
	      {
		// il faut déplacer le Bkp vers la gauche
		moveBkp=-1;
	      }                        
	  }
            
	if ((SmoothingNext <= NextLogRatio) && (NextLogRatio <= Smoothing))
	  {
	    if ( (NextLogRatio-SmoothingNext)>(Smoothing-NextLogRatio))
	      {
		// il faut déplacer le Bkp vers la droite
		moveBkp=1;
	      }
                
	  }

	if (LogRatio <= SmoothingNext)
	  {
	    moveBkp=-1;
	  }

	if (NextLogRatio>=Smoothing)
	  {
	    moveBkp=1;
	  }
      }
    // le créneau est plus bas à gauche qu'à droite
    else
      {
	if ((SmoothingNext >= LogRatio) && (LogRatio >= Smoothing))
	  {
	    if ((SmoothingNext-LogRatio) < (LogRatio - Smoothing))
	      {
		// il faut déplacer le Bkp vers la gauche
		moveBkp=-1;
	      }
	  }

	if ((SmoothingNext >= NextLogRatio) & (NextLogRatio >= Smoothing))
	  {
	    if ((SmoothingNext-NextLogRatio) > (NextLogRatio-Smoothing))
	      {
		// il faut déplacer le Bkp vers la droite
		moveBkp=1;
	      }
	  }

	if (LogRatio>=SmoothingNext)
	  {
	    moveBkp=-1;
	  }

	if (NextLogRatio<=Smoothing)
	  {
	    moveBkp=1;
	  }
      }
    return moveBkp;

  }






}
