/*****************************************************************************/
/* Copyright (C) 2004 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2004                            */
/* Contact: bioinfo-staff@curie.fr                                           */
/* It is strictly forbidden to transfer, use or re-use this code             */
/* or part of it without explicit written authorization from Institut Curie. */
/*****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <map>

#include "glad-function.h"
#include "glad.h"

#ifdef IS_MAC_OS
#include <limits.h>
#else
#include <values.h>
#endif

#ifndef MAXDOUBLE
#include <float.h>
#define MAXDOUBLE DBL_MAX
#endif

using namespace std;

extern "C"
{


  /*************************************/
  /* fonctions utilisées dans filterBkp */
  /*************************************/

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
		       const int *l)
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

  }

  void updateLevel (const int *Chromosome,
		    const int *Breakpoints,
		    int *Level,
		    //		    int *Region,
		    const int *PosOrder,
		    double *NextLogRatio,
		    const double *LogRatio,
		    const int *maxLevel,
		    const int *l)
  {
    const int nb=*l;
    int pos_moins_un;
    int pos;
    int idLevel=*maxLevel;
 
    for (pos=1;pos<nb;pos++)
      {
	pos_moins_un=pos-1;
	if (Chromosome[pos]==Chromosome[pos_moins_un])
	  {
	    if (Breakpoints[pos_moins_un]!=1)
	      {
		Level[pos]=Level[pos_moins_un];
		//		Region[pos]=Region[pos_moins_un];
	      }
	    if (Breakpoints[pos_moins_un]==1)
	      {
		NextLogRatio[pos_moins_un]=LogRatio[pos];
		if (Level[pos_moins_un]==Level[pos])
		  {
		    idLevel++;
		    Level[pos]=idLevel;
		  }
	      }
	  }
      }
  }


  void updateOutliers (int *OutliersAws,
		       int *Level,
		       int *Breakpoints,
		       double *Smoothing,
		       const int *l)
  {
  
    int pos;
    int pos_moins_un;
    const int nb=*l-1;

    for (pos=1;pos<nb;pos++)
      {
	pos_moins_un=pos-1;
	if (Level[pos_moins_un]==Level[pos+1] && Level[pos_moins_un]!=Level[pos])
	  {
	    Level[pos]=Level[pos_moins_un];
	    Breakpoints[pos_moins_un]=0;
	    Breakpoints[pos]=0;
	    OutliersAws[pos]=1;
	    Smoothing[pos]=Smoothing[pos_moins_un];
	  }
      }

  }

  /*************************************/
  /* fonctions utilisées dans MoveBkp */
  /*************************************/

  void MoveBkp_updateOutliers (int OutliersAws[],
			       int OutliersTot[],
			       int Level[],
			       int Region[],
			       int Breakpoints[],
			       double Smoothing[],
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


  void moveBkp (int *ZoneGNL,
		int *Level,
		int *Breakpoints,
		int *OutliersTot,
		int *OutliersAws,
		const int *Chromosome,
		int *RecomputeSmt,
		const int *l)
  {

    int pos;
    int pos_moins_un;
    int pos_plus_un;
    int pos_plus_deux;
    const int nb=*l;
    const int nb_moins_un=*l-1;
    const int nb_moins_deux=*l-2;

    for (pos=1;pos<nb;pos++)
      {
	pos_moins_un=pos-1;
	if (Chromosome[pos]==Chromosome[pos_moins_un])
	  {
	    pos_plus_un=pos+1;
	    if (OutliersTot[pos]!=0 && Breakpoints[pos]==1 && ZoneGNL[pos]==ZoneGNL[pos_plus_un] && ZoneGNL[pos_moins_un]!=ZoneGNL[pos_plus_un])
	      {
		*RecomputeSmt=1;
		Breakpoints[pos]=0;
		Breakpoints[pos_moins_un]=1;
		OutliersTot[pos]=0;
		OutliersAws[pos]=0;
		Level[pos]=Level[pos_plus_un];

	      }

	    if (pos < nb_moins_un && Breakpoints[pos]==1 && OutliersTot[pos_plus_un]!=0 && ZoneGNL[pos]==ZoneGNL[pos_plus_un] && ZoneGNL[pos_moins_un]!=ZoneGNL[pos+1])
	      {
		*RecomputeSmt=1;	  
		Breakpoints[pos]=0;
		Breakpoints[pos_plus_un]=1;
		OutliersTot[pos_plus_un]=0;
		OutliersAws[pos_plus_un]=0;
		Level[pos_plus_un]=Level[pos];

		if (pos < nb_moins_deux)
		  {
		    pos_plus_deux=pos+2;
		    if (Chromosome[pos_plus_un]==Chromosome[pos_plus_deux])
		      {
			if (Level[pos_plus_un]==Level[pos_plus_deux])
			  {
			    Breakpoints[pos_plus_un]=0;
			  }
		      }		    		        
		  }
	      }
	  }
      }
  }



  /*******************************************/
  /* fonctions utilisées dans chrBreakpoints */
  /*******************************************/

  void awsBkp (const double *Smoothing,
	       int *OutliersAws,
	       int *Level,
	       int *nbregion,
	       int *regionChr,
	       int *rupture,
	       int *bkp_detected,
	       const int *l)
  {

    int j;
    int j_moins_un;
    int j_plus_un;
    int last_bkp_pos = -1;
    const int nb = *l;
    const int nb_moins_un = *l - 1;

    for (j = 1; j < nb; j++)
      {              
	j_moins_un = j - 1;
	j_plus_un = j + 1;
	/*Outliers detection*/

	if (Smoothing[j_moins_un] != Smoothing[j] && 
	    Smoothing[j] != Smoothing[j_plus_un] &&
	    Smoothing[j_moins_un] == Smoothing[j_plus_un] &&
	    j < nb_moins_un)
	  {
	    if (OutliersAws[j_moins_un] == 0)
	      {
		if (Smoothing[j] > Smoothing[j_plus_un])
		  {
		    OutliersAws[j] = 1;
		    Level[j] = Level[j_moins_un];
		  }

		else
		  {
		    OutliersAws[j] = -1;
		    Level[j] = Level[j_moins_un];
		  }
	      }
                  
	    regionChr[j] = *nbregion; /*en entrée, on a region[j==0]=*nbregion*/
	  }
              
	/* not Outliers*/
	/* Pour les outliers, leur Level correspond à celui de la région */
	/* de laquelle ils sont issus. Cela permet de faciliter l'affectation */
	/* de leur statut dans la fonction affectationGNL  */ 
            
	else
	  {
	    if (Smoothing[j_moins_un] != Smoothing[j] &&
		OutliersAws[j_moins_un] == 0)
	      {
		if (j == 1 || j == nb_moins_un)
		  {
		    regionChr[j] = *nbregion;

		    if (j == 1)
		      {
			if (Smoothing[0] > Smoothing[1])
			  {
			    OutliersAws[0] = 1;
			    Level[0] = Level[1];
			  }

			else
			  {
			    OutliersAws[0] = -1;
			    Level[0] = Level[1];
			  }
		      }
                          
		    else
		      {
			if (Smoothing[j] > Smoothing[j_moins_un])
			  {                                  
			    OutliersAws[j] = 1;
			    Level[j] = Level[j_moins_un];
			  }

			else
			  {
			    OutliersAws[j] = -1;
			    Level[j] = Level[j_moins_un];
			  }
		      }
		  }
                      
		else
		  {
		    if (last_bkp_pos == j_moins_un)
		      {
			if (Smoothing[j_moins_un] > Smoothing[j_moins_un -1])
			  {                                  
			    OutliersAws[j_moins_un] = 1;
			  }
			else
			  {
			    OutliersAws[j_moins_un] = -1;
			  }

			// on supprime le Bkp précédent
			rupture[j_moins_un] = 0;

			Level[j_moins_un] = Level[j_moins_un - 1];
			regionChr[j_moins_un] = regionChr[j_moins_un - 1];


			// on ajoute un Bkp
			// il ne faut pas incrémenter *nbregion car cela aura déjà été fait
			regionChr[j] = *nbregion;
			rupture[j] = 1;
			last_bkp_pos=j;
			*bkp_detected = 1;
		      }
		    else
		      {
			*nbregion = *nbregion + 1;
			regionChr[j] = *nbregion;
			rupture[j] = 1;
			last_bkp_pos=j;
			*bkp_detected = 1;
		      }
		  }
	      }

	    else
	      {
		regionChr[j] = *nbregion;
	      }
	  }
      }
  }

  

  void updateBkpRL (int *Region,
		    int *OutliersAws,
		    int *Breakpoints,
		    const int *Chromosome,
		    const int *PosOrder,
		    double *NextLogRatio,
		    const double *LogRatio,
		    const int *l)
  {

    int i;
    int i_moins_un;
    int i_plus_un;
    int i_moins_deux;
    const int nb=*l;
    const int nb_moins_un=*l-1;
    const int nb_moins_deux=*l-2;

    for (i=1;i<nb;i++)
      {
	i_moins_un=i-1;
	if (i==1 || i==nb_moins_un)
	  {
	    if (Region[i]!=Region[i_moins_un])
	      {
		if (i==1)
		  {
		    OutliersAws[0]=1;
		    Region[0]=Region[1];
		  }
		else
		  {
		    OutliersAws[nb_moins_un]=1;
		    Region[nb_moins_un]=Region[nb_moins_deux];
		  }
	      }
	  }
	else
	  {
	    i_plus_un=i+1;
	    if (Chromosome[i]!=Chromosome[i_moins_un])
	      {
		i_moins_deux=i-2;
		if (Region[i_moins_un]!=Region[i_moins_deux])
		  {
		    Region[i_moins_un]=Region[i_moins_deux];
		    OutliersAws[i_moins_un]=1;
		  }

		if (Region[i_plus_un]!=Region[i])
		  {
		    Region[i_plus_un]=Region[i];
		    OutliersAws[i]=1;
		  }                
	      }
	    else
	      {
		if (Region[i]!=Region[i_moins_un] && Region[i_plus_un]!=Region[i] && Region[i_plus_un]==Region[i_moins_un])
		  {
		    if (OutliersAws[i_moins_un]==0)
		      {
			OutliersAws[i]=1;
			Region[i]=Region[i_moins_un];
		      }
		  }
		else
		  {
		    if (Region[i]!=Region[i_moins_un] && OutliersAws[i_moins_un]==0)
		      {
			Breakpoints[i_moins_un]=1;
			NextLogRatio[i_moins_un]=LogRatio[i];
		      }
		  }                                               
	      }
	  }
      }
  }



  void MoveBkp_Delete_Bkp(const int subBkpInfo_MoveBkp[],
			  const int subBkpInfo_PosOrder[],
			  int Breakpoints[],
			  int OutliersTot[],
			  int OutliersAws[],
			  int OutliersMad[],
			  int Level[],
			  int Region[],
			  double Smoothing[],
			  int GNL[],
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
    const int nb=*l;

    *minG=MAXDOUBLE;
    *minAmp=MAXDOUBLE;
    *maxL=-MAXDOUBLE;
    *maxDel=-MAXDOUBLE;

    for (i=0;i<nb;i++)
      {
	if (OutliersTot[i]==0)
	  {

	    switch(ZoneGNL[i])
	      {
	      case 0:
		break;

	      case 1:
		if(Smoothing[i] < *minG)
		  {
		    *minG=Smoothing[i];
		  }
		break;

	      case -1:
		if (Smoothing[i] > *maxL)
		  { 
		    *maxL=Smoothing[i]; 
		  }
		break;

	      case 2:
		if(Smoothing[i] < *minAmp)
		  {
		    *minAmp=Smoothing[i];
		  }
		break;

	      case -10:
		if(Smoothing[i] > *maxDel)
		  {
		    *maxDel=Smoothing[i];
		  }
		break;

	      }
	  }

      }
  }


  void MoveBkp_updateGNL(int ZoneGNL[],
			 const double Smoothing[],
			 const int OutliersTot[],
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


  /*************************************/
  /* fonctions utilisées dans daglad   */
  /*************************************/

  void delete_contiguous_bkp(int *BkpInfo_BkpToDel,
			     double *BkpInfo_Gap,
			     double *BkpInfo_LogRatio,
			     int *BkpInfo_NextPosOrder,
			     int *BkpInfo_PosOrder,
			     int *BkpInfo_Side,
			     double *BkpInfo_Sigma,
			     double *BkpInfo_Smoothing,
			     double *BkpInfo_SmoothingNext,
			     double *BkpInfo_Weight,
			     int *nb_Bkp,
			     int *RecomputeGNL,
			     const int *nbsigma)
  {

    // nb_Bkp=length(profileCGH$BkpInfo[,1])
    int i;
    int i_moins_un;
    const int l=*nb_Bkp;
    double UpLeft, LowLeft, UpRight, LowRight, DiffLeft, DiffRight;
    double LRV;
    double SigChr;


    //  for (i in 2:length(profileCGH$BkpInfo[,1]))
    for (i=1;i<l;i++)
      {
	i_moins_un=i-1;
	if (BkpInfo_PosOrder[i]==BkpInfo_NextPosOrder[i_moins_un] && BkpInfo_BkpToDel[i_moins_un]==0)
	  {
	    SigChr=BkpInfo_Sigma[i];
	    // on regarde d'abord à gauche
	    UpLeft=BkpInfo_Smoothing[i_moins_un] + 3*SigChr;
	    LowLeft=BkpInfo_Smoothing[i_moins_un] - 3*SigChr;

	    // On regarde ce qui se passe à droite
	    UpRight=BkpInfo_SmoothingNext[i] + 3*SigChr;
	    LowRight=BkpInfo_SmoothingNext[i] - 3*SigChr;

	    LRV=BkpInfo_LogRatio[i];
                    
	    if (((LRV > LowLeft) && (LRV < UpLeft)) || ((LRV > LowRight) && (LRV < UpRight)))
	      {
		*RecomputeGNL=1;
		if (((LRV > LowLeft) && (LRV < UpLeft)) && ((LRV > LowRight) && (LRV < UpRight)))
		  {
		    // attention, lors de la suppression d'un Bkp, le gap n'est plus bon
		    // d'où recalcul du weight
		    // je ne vois pas où il est recalculé!!!!
		    DiffLeft=fabs(LRV - BkpInfo_Smoothing[i_moins_un]);
		    DiffRight=fabs(LRV - BkpInfo_SmoothingNext[i]);
		    if (DiffRight<DiffLeft)
		      {
			// On fusionne le Bkp avec la région à droite
			BkpInfo_BkpToDel[i]=1;
			BkpInfo_Side[i]=1;
			BkpInfo_Gap[i_moins_un]=fabs(BkpInfo_Smoothing[i_moins_un]-BkpInfo_SmoothingNext[i]);
			BkpInfo_Weight[i_moins_un]=1 - kernelpen(BkpInfo_Gap[i_moins_un], *nbsigma*BkpInfo_Sigma[i_moins_un]);
                                
		      }
		    else
		      {
			// On fusionne le Bkp avec la région à gauche
			BkpInfo_BkpToDel[i_moins_un]=1;
			BkpInfo_Side[i_moins_un]=0;
			BkpInfo_Gap[i]=fabs(BkpInfo_Smoothing[i_moins_un]-BkpInfo_SmoothingNext[i]);
			BkpInfo_Weight[i]=1 - kernelpen(BkpInfo_Gap[i_moins_un], *nbsigma*BkpInfo_Sigma[i_moins_un]);
                                

		      }
		  }
		else
		  {
		    if (((LRV > LowLeft) && (LRV < UpLeft)))
		      {
			// On fusionne le Bkp avec la région à gauche
			BkpInfo_BkpToDel[i_moins_un]=1;
			BkpInfo_Side[i_moins_un]=0;
			BkpInfo_Gap[i]=fabs(BkpInfo_Smoothing[i_moins_un]-BkpInfo_SmoothingNext[i]);
			BkpInfo_Weight[i]=1 - kernelpen(BkpInfo_Gap[i_moins_un], *nbsigma*BkpInfo_Sigma[i_moins_un]);
                                
		      }
		    else
		      {
			// On fusionne le Bkp avec la région à droite
			BkpInfo_BkpToDel[i]=1;
			BkpInfo_Side[i]=1;
			BkpInfo_Gap[i_moins_un]=fabs(BkpInfo_Smoothing[i_moins_un]-BkpInfo_SmoothingNext[i]);
			BkpInfo_Weight[i_moins_un]=1 - kernelpen(BkpInfo_Gap[i_moins_un], *nbsigma*BkpInfo_Sigma[i_moins_un]);                            
		      }
		  }                                        
	      }
	  }
      }
  }




  /*************************************/
  /* fonctions utilisées dans BkpInfo.R  */
  /*************************************/

  void make_BkpInfo(const double *BkpInfo_Gap,
		    int *BkpInfo_GNLchange,
		    const double *BkpInfo_Value,
		    double *BkpInfo_Weight,
		    int *BkpInfo_ZoneGNL,
		    const int *BkpInfo_ZoneGNLnext,
		    const int *nb_Bkp,
		    const double *nbsigma)
  {

    const int l=*nb_Bkp;
    int i;

    for (i=0;i<l;i++)
      {
	BkpInfo_Weight[i]=1 - kernelpen(BkpInfo_Gap[i], *nbsigma*BkpInfo_Value[i]);
	if (BkpInfo_ZoneGNL[i]==BkpInfo_ZoneGNLnext[i])
	  {
	    BkpInfo_GNLchange[i]=0;
	  }
	else
	  {
	    BkpInfo_GNLchange[i]=1;
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
    return(moveBkp);

  }


  /*******************************************/
  /* fonctions utilisées dans OutliersGNL.R  */
  /*******************************************/

  void OutliersGNL(int * OutliersTot,
		   int *ZoneGNL,
		   const double *LogRatio,
		   const double * Smoothing,
		   const double *seuilsupValue,
		   const double *seuilinfValue,
		   const double *ampliconValue,
		   const double *deletionValue,
		   const double *NormalRefValue,
		   const int *l)
  {
    int i;
    const int nb=*l;

    const double seuilsup=*seuilsupValue;
    const double seuilinf=*seuilinfValue;
    const double amplicon=*ampliconValue;
    const double deletion=*deletionValue;
    const double NormalRef=*NormalRefValue;

    int checkGain=0;
    int checkLost=0;
    int checkNormal=0;
    int checkAlert=0;

    double LogRatio_moins_NormalRef;
    double minNormal=MAXDOUBLE;
    double maxNormal=-MAXDOUBLE;
    double minGain=MAXDOUBLE;
    double maxLost=-MAXDOUBLE;

    /*   printf("seuilsup=%f\n",seuilsup); */
    /*   printf("seuilinf=%f\n",seuilinf); */
    /*   printf("amplicon=%f\n",amplicon); */
    /*   printf("deletion=%f\n",deletion); */

    for (i=0;i<nb;i++)
      {
	//////////////////////////
	// On regarde les Outliers
	//////////////////////////
	if(OutliersTot[i]!=0)
	  {
	    // On met le GNL de tous les Outliers à 0
	    ZoneGNL[i]=0;

	    // Calcul de la différence entre le LogRatio et NormalRef
	    if(NormalRef!=0)
	      {
		LogRatio_moins_NormalRef=LogRatio[i]-NormalRef;
	      }
	    else
	      {
		LogRatio_moins_NormalRef=LogRatio[i];
	      }

	    // Gain et Amplicon
	    if(LogRatio_moins_NormalRef>seuilsup)
	      {
		// On a un Amplicon
		if(LogRatio_moins_NormalRef >= amplicon)
		  {
		    ZoneGNL[i]=2;
		  }
		// On a un Gain
		else
		  {
		    ZoneGNL[i]=1;
		  }
	      }
	    // Perte et Deletion
	    else
	      {
		if(LogRatio_moins_NormalRef<seuilinf)
		  {
		    // On a une deletion
		    if(LogRatio_moins_NormalRef<deletion)
		      {
			ZoneGNL[i]=-10;
		      }
		    // On a une Perte
		    else
		      {
			ZoneGNL[i]=-1;
		      }
		  }
	      }
	  }
	/////////////////////////////////////////////////////
	// Récupération des min/max pour le Normal/Gain/Perte
	/////////////////////////////////////////////////////
	else
	  {
	    switch(ZoneGNL[i])
	      {
	      case 0:
		if(Smoothing[i]<minNormal)
		  {
		    minNormal=Smoothing[i];
		  }
		if(Smoothing[i]>maxNormal)
		  {
		    maxNormal=Smoothing[i];
		  }
		checkNormal=1;
		break;

	      case 1:
		if(Smoothing[i]<minGain)
		  {
		    minGain=Smoothing[i];
		  }
		checkGain=1;
		break;

	      case -1:
		if(Smoothing[i]>maxLost)
		  {
		    maxLost=Smoothing[i];
		  }
		checkLost=1;
		break;
	      }
	  }
      }


    ////////////////////////////////////////////////////////////////////
    // On fait une seconde boucle pour mettre à jour le GNL des Outliers 
    // Et vérifier la cohérence des valeurs
    ////////////////////////////////////////////////////////////////////

    for(i=0;i<nb;i++)
      {
	//////////////////////////
	// On regarde les Outliers
	//////////////////////////
	if(OutliersTot[i]!=0)
	  {
	    if(ZoneGNL[i]==0)
	      {
		// comparaison avec minGain et maxLost
		if(LogRatio[i]>minGain)
		  {
		    ZoneGNL[i]=1;
		  }
		else
		  {
		    if(LogRatio[i]<maxLost)
		      {
			ZoneGNL[i]=-1;
		      }
		  }
	      }
	  }
	else
	  {
	    if(checkLost && checkGain)
	      {
		if(checkNormal)
		  {
		    if (maxLost>minNormal)
		      {
			if(Smoothing[i]<=maxLost && ZoneGNL[i]==0)
			  {
			    ZoneGNL[i]=-1;
			    checkAlert=1;
			  }
		      }
		    if(minGain<minNormal)
		      {
			if(Smoothing[i]>=minGain && ZoneGNL[i]==0)
			  {
			    ZoneGNL[i]=1;
			    checkAlert=1;
			  }
		      }
		  }
	      }
	  }
      }

    if (maxLost>minGain)
      {
	checkAlert=1;
      }

    if (checkAlert)
      {
	printf("In function OutliersGNL: Inconsistency for smoothing values vs. GNL status has been corrected)\n");
      }
  }




  /********************************************************************/
  /* fonction pour faire une aggregation avec le calcul de la médian  */
  /********************************************************************/


  void compute_median_smoothing(const double LogRatio[],
				const int ByValue[],
				double Smoothing[],
				const int *l)
  {
    int i,j;
    const int nb=*l;
    double *median_ByValue;
    int *unique_ByValue;
    int nb_unique_ByValue;

    map<int, vector<double> > agg_LogRatio;
    map<int, vector<double> >::iterator it_agg_LogRatio;

    for(i = 0; i < nb; i++)
      {
	agg_LogRatio[ByValue[i]].push_back(LogRatio[i]);
      }


    median_ByValue = (double *)malloc(agg_LogRatio.size() * sizeof(double));
    unique_ByValue = (int *)malloc(agg_LogRatio.size() * sizeof(int));
    it_agg_LogRatio=agg_LogRatio.begin();

    for (j = 0; j < agg_LogRatio.size(); j++)
      {
	median_ByValue[j] = median_vector_double(it_agg_LogRatio->second);
	unique_ByValue[j] = it_agg_LogRatio->first;

	it_agg_LogRatio++;

      }

    nb_unique_ByValue = (int)agg_LogRatio.size();
    my_merge(ByValue,
	     Smoothing,
	     unique_ByValue,
	     median_ByValue,
	     l,
	     &nb_unique_ByValue);




    free(median_ByValue);
    free(unique_ByValue);
  }

  /**********************************************************/
  /* fonction pour le calcul des clusters Loss/Normal/Gain  */
  /**********************************************************/


  void compute_cluster_LossNormalGain(// variables pour faire la jointure
				      const int ZoneGen[],
				      int value_dest[],
				      const int *length_dest,
				      const double Smoothing[],
				      const double *forceGL1Value,
				      const double *forceGL2Value,
				      const double *NormalRefValue,
				      const double *ampliconValue,
				      const double *deletionValue,
				      //variables pour calcul la médiane par cluster
				      const double LogRatio[],
				      const int NormalRange[])
  {

    int i,j;
    int nb = *length_dest;
    int NormalCluster;
    int NormalCluster_not_detected = 1;

    int *MedianCluster_ZoneGen;
    int *MedianCluster_ZoneGNL;
    double *MedianCluster_Median;
    int nb_unique_ZoneGen;

    double RefNorm;
    vector<int>::iterator it_new_end_NormalCluster;


    map<int, vector<double> > agg_LogRatio;
    map<int, vector<double> >::iterator it_agg_LogRatio;


    // On récupére les valeurs de LogRatio pour chaque ZoneGen
    for (i = 0; i < nb; i++)
      {
	agg_LogRatio[ZoneGen[i]].push_back(LogRatio[i]);

	// le cluster correspondant au normal est celui qui comprend
	// le NormalRange 0
	if(NormalRange[i] == 0 && NormalCluster_not_detected)
	  {
	    NormalCluster = ZoneGen[i];
	    NormalCluster_not_detected = 0;
	  }
      }

    // On calcule la médiane par ZoneGen
    MedianCluster_Median = (double *)malloc(agg_LogRatio.size() * sizeof(double));
    MedianCluster_ZoneGen = (int *)malloc(agg_LogRatio.size() * sizeof(int));
    MedianCluster_ZoneGNL = (int *)malloc(agg_LogRatio.size() * sizeof(int));
    it_agg_LogRatio=agg_LogRatio.begin();

    for (j = 0; j < agg_LogRatio.size(); j++)
      {
	MedianCluster_Median[j] = median_vector_double(it_agg_LogRatio->second);
	MedianCluster_ZoneGen[j] = it_agg_LogRatio->first;

	if(NormalCluster == it_agg_LogRatio->first)
	  {
	    RefNorm = MedianCluster_Median[j];
	  }

	it_agg_LogRatio++;

      }

    for (j = 0; j < agg_LogRatio.size(); j++)
      {
	MedianCluster_ZoneGNL[j] = 0;

	if(MedianCluster_Median[j] > RefNorm)
	  {
	    MedianCluster_ZoneGNL[j] = 1;
	  }
	else
	  {
	    if(MedianCluster_Median[j] < RefNorm)
	      {
		MedianCluster_ZoneGNL[j] = -1;
	      }
	  }
      }

    nb_unique_ZoneGen = (int)agg_LogRatio.size();
    my_merge_int_forceGL(ZoneGen,
			 value_dest,
			 MedianCluster_ZoneGen,
			 MedianCluster_ZoneGNL,
			 length_dest,
			 &nb_unique_ZoneGen,
			 Smoothing,
			 forceGL1Value,
			 forceGL2Value,
			 NormalRefValue,
			 ampliconValue,
			 deletionValue);


    free(MedianCluster_ZoneGen);
    free(MedianCluster_Median);
    free(MedianCluster_ZoneGNL);

  }

  void compute_NormalRange(const double Smoothing[],
			   const double *NormalRef,
			   const int Level[],
			   int NormalRange[],
			   const double *deltaN,
			   const int *l)
  {

    int i;
    const int nb = *l;

    for (i = 0; i < nb; i++)
      {
	if(fabs(Smoothing[i] - *NormalRef) <= *deltaN)
	  {
	    NormalRange[i] = 0;
	  }
	else
	  {
	    NormalRange[i] = Level[i];
	  }
      }

  }

}





