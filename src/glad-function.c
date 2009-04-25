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

#include "glad-function.h"

#ifdef IS_MAC_OS
#include <limits.h>
#else
#include <values.h>
#endif

#ifndef MAXDOUBLE
#include <float.h>
#define MAXDOUBLE DBL_MAX
#endif






/*************************************/
/* fonctions utilisées dans filterBkp */
/*************************************/

void updateLevel (int *Chromosome,
                  int *Breakpoints,
                  int *Level,
                  int *PosOrder,
		  double *NextLogRatio,
		  double *LogRatio,
		  int *maxLevel,
                  int *l)
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
                     int *l)
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

void updateOutliersMoveBkp (int *OutliersAws,
			    int *OutliersTot,
			    int *Level,
			    int *Breakpoints,
			    double *Smoothing,
			    int *ZoneGNL,
			    int *l)
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
	  ZoneGNL[pos]=ZoneGNL[pos_moins_un];
	  Smoothing[pos]=Smoothing[pos_moins_un];
	}

    }

}


void moveBkp (int *ZoneGNL,
              int *Level,
              int *Breakpoints,
              int *OutliersTot,
              int *OutliersAws,
	      int *Chromosome,
              int *RecomputeSmt,
              int *l)
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

void awsBkp (double *Smoothing,
             int *OutliersAws,
             int *Level,
             int *nbregion,
             int *regionChr,
             int *rupture,
             int *l)
{

  int j;
  int j_moins_un;
  int j_plus_un;
  const int nb=*l;
  const int nb_moins_un=*l-1;

  for (j=1;j<nb;j++)
    {              
      j_moins_un=j-1;
      j_plus_un=j+1;
      /*Outliers detection*/

      if (Smoothing[j]!=Smoothing[j_moins_un] && Smoothing[j_plus_un]!=Smoothing[j] && Smoothing[j_plus_un]==Smoothing[j_moins_un] && j<nb_moins_un)
	{
	  if (OutliersAws[j_moins_un]==0)
	    {
	      if (Smoothing[j]>Smoothing[j_plus_un])
		{
		  OutliersAws[j]=1;
		  Level[j]=Level[j_moins_un];
		}

	      else
		{
		  OutliersAws[j]=-1;
		  Level[j]=Level[j_moins_un];
		}
	    }
                  
	  regionChr[j]=*nbregion; /*en entrée, on a region[j==0]=*nbregion*/
	  rupture[j]=0;
	}
              
      /* not Outliers*/
      /* Pour les outliers, leur Level correspond à celui de la région */
      /* de laquelle ils sont issus. Cela permet de faciliter l'affectation */
      /* de leur statut dans la fonction affectationGNL  */ 
            
      else
	{
	  if (Smoothing[j]!=Smoothing[j_moins_un] && OutliersAws[j_moins_un]==0)
	    {
	      if (j==1 || j==nb_moins_un)
		{
		  regionChr[j]=*nbregion;
		  rupture[j]=0;

		  if (j==1)
		    {
		      if (Smoothing[0]>Smoothing[1])
			{
			  OutliersAws[0]=1;
			  Level[0]=Level[1];
			}

		      else
			{
			  OutliersAws[0]=-1;
			  Level[0]=Level[1];
			}
		    }
                          
		  else
		    {
		      if (Smoothing[j]>Smoothing[j_moins_un])
			{                                  
			  OutliersAws[j]=1;
			  Level[j]=Level[j_moins_un];
			}

		      else
			{
			  OutliersAws[j]=-1;
			  Level[j]=Level[j_moins_un];
			}
		    }
		}
                      
	      else
		{
		  *nbregion=*nbregion+1;
		  regionChr[j]=*nbregion;
		  rupture[j]=1;
		}
	    }

	  else
	    {
	      regionChr[j]=*nbregion;
	      rupture[j]=0;
	    }
	}

    }


}

  

void updateBkpRL (int *Region,
                  int *OutliersAws,
                  int *Breakpoints,
                  int *Chromosome,
                  int *PosOrder,
		  double *NextLogRatio,
		  double *LogRatio,
                  int *l)
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



void loopMoveBkp(int *subBkpInfo_MoveBkp,
		 int *subBkpInfo_PosOrder,
		 int *CGH_Breakpoints,
		 int *CGH_OutliersTot,
		 int *CGH_OutliersAws,
		 int *CGH_OutliersMad,
		 int *CGH_Level,
		 double *CGH_Smoothing,
		 int *CGH_GNL,
		 int *l)
{

  int i;
  int indexPos;
  int indexPosNext;
  int indexPosBefore;
  const int nb=*l;

  for (i=0;i<nb;i++)
    {


      if (subBkpInfo_MoveBkp[i]==1)
	{
	  // déplacement à droite
	  indexPos=subBkpInfo_PosOrder[i]-1;
	  indexPosNext=indexPos+1;
	  CGH_Breakpoints[indexPos]=0;
	  CGH_Breakpoints[indexPosNext]=1;
	  CGH_OutliersTot[indexPosNext]=0;
	  CGH_OutliersAws[indexPosNext]=0;
	  CGH_OutliersMad[indexPosNext]=0;
	  CGH_OutliersTot[indexPos]=0;
	  CGH_OutliersAws[indexPos]=0;
	  CGH_OutliersMad[indexPos]=0;
	  CGH_Level[indexPosNext]=CGH_Level[indexPos];
	  CGH_Smoothing[indexPosNext]=CGH_Smoothing[indexPos];
	  CGH_GNL[indexPosNext]=CGH_GNL[indexPos];
   
	}
      if (subBkpInfo_MoveBkp[i]==-1)
	{
	  // déplacement à gauche
	  indexPos=subBkpInfo_PosOrder[i]-1;
	  indexPosBefore=indexPos-1;
	  indexPosNext=indexPos+1;
	  CGH_Breakpoints[indexPos]=0;
	  CGH_Breakpoints[indexPosBefore]=1;
	  CGH_OutliersTot[indexPosBefore]=0;
	  CGH_OutliersAws[indexPosBefore]=0;
	  CGH_OutliersMad[indexPosBefore]=0;
	  CGH_OutliersTot[indexPos]=0;
	  CGH_OutliersAws[indexPos]=0;
	  CGH_OutliersMad[indexPos]=0;
	  CGH_Level[indexPos]=CGH_Level[indexPosNext];
	  CGH_Smoothing[indexPos]=CGH_Smoothing[indexPosNext];
	  CGH_GNL[indexPos]=CGH_GNL[indexPosNext];
                    
	}
                
    }
}




void rangeGainLoss(double *Smoothing,
		   int *ZoneGNL,
		   int *OutliersTot,
		   double *minG,
		   double *maxL,
		   double *minAmp,
		   double *maxDel,
		   int *l)
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


void updateGNL(int *ZoneGNL,
	       double *Smoothing,
	       int *OutliersTot,
	       int *l)
{
  double *minG;
  double *maxL;
  double *minAmp;
  double *maxDel;
  int i;
  const int nb=*l;

  minG=malloc(1 * sizeof(double));
  maxL=malloc(1 * sizeof(double));
  minAmp=malloc(1 * sizeof(double));
  maxDel=malloc(1 * sizeof(double));

  rangeGainLoss(Smoothing, ZoneGNL, OutliersTot, minG, maxL, minAmp, maxDel, l);


  for (i=0; i<nb;i++)
    {
      ZoneGNL[i]=0;

      if (Smoothing[i]>=*minG || Smoothing[i]>=*minAmp)
	{
	  ZoneGNL[i]=1;
	 
	  if(Smoothing[i]>=*minAmp)
	    {
	      ZoneGNL[i]=2;
	    }

	}

      else
	{
	  if (Smoothing[i]<=*maxL || Smoothing[i]<=*maxDel)
	    {
	      ZoneGNL[i]=-1;

	      if(Smoothing[i]<=*maxDel)
		{
		  ZoneGNL[i]=-10;
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
			   int *nbsigma)
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
/* fonctions utilisées dans Kernel.R  */
/*************************************/


double kernelpen(double value, const double d)
{
  double tricubic;
  if (value>d)
    return(0);

  // index <- which(x<=param["d"])
  // k[index] <- (1-(x[index]/param["d"])^3)^3

  tricubic=value/d;
  tricubic=tricubic*tricubic*tricubic;
  tricubic=1-tricubic;

  return(tricubic*tricubic*tricubic);

}


/*************************************/
/* fonctions utilisées dans BkpInfo.R  */
/*************************************/

void make_BkpInfo(double *BkpInfo_Gap,
		  int *BkpInfo_GNLchange,
		  double *BkpInfo_Value,
		  double *BkpInfo_Weight,
		  int *BkpInfo_ZoneGNL,
		  int *BkpInfo_ZoneGNLnext,
		  int *nb_Bkp,
		  double *nbsigma)
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
	      if(LogRatio_moins_NormalRef>amplicon)
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


/*   printf("minNormal=%f\n",minNormal); */
/*   printf("maxNormal=%f\n",maxNormal); */
/*   printf("minGain=%f\n",minGain); */
/*   printf("maxLost=%f\n",maxLost); */

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


/*********************/
/* fonctions de tri  */
/*********************/

void my_order_int(int* array, int *order, const int *leftValue, const int *rightValue)
{
  int left=*leftValue;
  int right=*rightValue;

  quicksort_int(array, order, left, right);

}

//Quicksort the array
void quicksort_int(int* array, int *order, int left, const int right)
{
  if(left >= right)
    return;
 
  int index = partition(array, order, left, right);
  quicksort_int(array, order, left, index - 1);
  quicksort_int(array, order, index + 1, right);
}
 
//Partition the array into two halves and return the
//index about which the array is partitioned
int partition(int* array, int *order, int left, int right)
{
  //Makes the leftmost element a good pivot,
  //specifically the median of medians
  findMedianOfMedians(array, order, left, right);
  int pivotIndex = left, pivotValue = array[pivotIndex], index = left, i;
 
  swap(array, order, pivotIndex, right);
  for(i = left; i < right; i++)
    {
      if(array[i] < pivotValue)
        {
	  swap(array, order, i, index);
	  index += 1;
        }
    }
  swap(array, order, right, index);
 
  return index;
}
 
//Computes the median of each group of 5 elements and stores
//it as the first element of the group. Recursively does this
//till there is only one group and hence only one Median
int findMedianOfMedians(int* array, int *order, int left, int right)
{
  if(left == right)
    return array[left];
 
  int i, shift = 1;
  while(shift <= (right - left))
    {
      for(i = left; i <= right; i+=shift*5)
        {
	  int endIndex = (i + shift*5 - 1 < right) ? i + shift*5 - 1 : right;
	  int medianIndex = findMedianIndex(array, order, i, endIndex, shift);
 
	  swap(array, order, i, medianIndex);
        }
      shift *= 5;
    }
 
  return array[left];
}
 
//Find the index of the Median of the elements
//of array that occur at every "shift" positions.
int findMedianIndex(int* array, int *order, int left, int right, int shift)
{
  int i, groups = (right - left)/shift + 1, k = left + groups/2*shift;
  for(i = left; i <= k; i+= shift)
    {
      int minIndex = i, minValue = array[minIndex], j;
      for(j = i; j <= right; j+=shift)
	if(array[j] < minValue)
	  {
	    minIndex = j;
	    minValue = array[minIndex];
	  }
      swap(array, order, i, minIndex);
    }
 
  return k;
}
 
//Swap integer values by array indexes
void swap(int *array, int *order, int a, int b)
{
  int tmp  = array[a];
  int tmp_order=order[a];

  array[a] = array[b];
  array[b] = tmp;

  order[a] = order[b];
  order[b] = tmp_order;
}

