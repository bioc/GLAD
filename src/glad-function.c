/*****************************************************************************/
/*                                                                           */
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


/* void makeRegion (int *Level, */
/* 		 int *Chromosome, */
/* 		 int *ResLevel, */
/* 		 int *l) */
/* { */
/*   int pos; */
/*   int idLevel=1; */
 
/*   ResLevel[0]=idLevel; */
 
/*   for (pos=1;pos<*l;pos++) */
/*     { */
/*       if (Chromosome[pos]==Chromosome[pos-1]) */
/* 	{ */
/* 	  if (Level[pos]==Level[pos-1]) */
/* 	    { */
/* 	      ResLevel[pos]=ResLevel[pos-1]; */
/* 	    } */
/* 	  else */
/* 	    { */
/* 	      idLevel++; */
/* 	      ResLevel[pos]=idLevel; */
/* 	    } */
/* 	} */
/*       else */
/* 	{ */
/* 	  idLevel++; */
/* 	  ResLevel[pos]=idLevel; */
/* 	} */
/*     } */

/* } */





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
  int pos;
  int idLevel=*maxLevel;
 
/*   Level[0]=idLevel; */
 
  for (pos=1;pos<*l;pos++)
    {
      if (Chromosome[pos]==Chromosome[pos-1])
	{
	  if (Breakpoints[pos-1]!=1)
	    {
	      Level[pos]=Level[pos-1];
	    }
	  if (Breakpoints[pos-1]==1)
	    {
	      NextLogRatio[pos-1]=LogRatio[pos];
	      if (Level[pos-1]==Level[pos])
		{
		  idLevel++;
		  Level[pos]=idLevel;
		}

	    }
	}
/*       else */
/* 	{ */
/* 	  idLevel++; */
/* 	  Level[pos]=idLevel; */
/* 	} */
    }

}


void updateOutliers (int *OutliersAws,
                     int *Level,
                     int *Breakpoints,
		     double *Smoothing,
                     int *l)
{
  
  int pos;
  for (pos=1;pos<(*l-1);pos++)
    {

      if (Level[pos-1]==Level[pos+1] && Level[pos-1]!=Level[pos])
	{
	  Level[pos]=Level[pos-1];
	  Breakpoints[pos-1]=0;
	  Breakpoints[pos]=0;
	  OutliersAws[pos]=1;
	  Smoothing[pos]=Smoothing[pos-1];
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
  for (pos=1;pos<(*l-1);pos++)
    {
      
      if (Level[pos-1]==Level[pos+1] && Level[pos-1]!=Level[pos])
	{
	  Level[pos]=Level[pos-1];
	  Breakpoints[pos-1]=0;
	  Breakpoints[pos]=0;
	  OutliersAws[pos]=1;
	  ZoneGNL[pos]=ZoneGNL[pos-1];
	  Smoothing[pos]=Smoothing[pos-1];
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

  for (pos=1;pos<*l;pos++)
    {
      if (Chromosome[pos]==Chromosome[pos-1])
	{
	  if (OutliersTot[pos]!=0 && Breakpoints[pos]==1 && ZoneGNL[pos]==ZoneGNL[pos+1] && ZoneGNL[pos-1]!=ZoneGNL[pos+1])
	    {
	      *RecomputeSmt=1;
	      Breakpoints[pos]=0;
	      Breakpoints[pos-1]=1;
	      OutliersTot[pos]=0;
	      OutliersAws[pos]=0;
	      Level[pos]=Level[pos+1];

	    }

	  if (pos < (*l-1) && Breakpoints[pos]==1 && OutliersTot[pos+1]!=0 && ZoneGNL[pos]==ZoneGNL[pos+1] && ZoneGNL[pos-1]!=ZoneGNL[pos+1])
	    {
	      *RecomputeSmt=1;	  
	      Breakpoints[pos]=0;
	      Breakpoints[pos+1]=1;
	      OutliersTot[pos+1]=0;
	      OutliersAws[pos+1]=0;
	      Level[pos+1]=Level[pos];

	      if (pos < (*l-2))
		{

		  if (Chromosome[pos+1]==Chromosome[pos+2])
		    {
		      if (Level[pos+1]==Level[pos+2])
			{
			  Breakpoints[pos+1]=0;
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
  for (j=1;j<*l;j++)
    {
              
      /*Outliers detection*/

      if (Smoothing[j]!=Smoothing[j-1] && Smoothing[j+1]!=Smoothing[j] && Smoothing[j+1]==Smoothing[j-1] && j<(*l-1))
	{
	  if (OutliersAws[j-1]==0)
	    {
	      if (Smoothing[j]>Smoothing[j+1])
		{
		  OutliersAws[j]=1;
		  Level[j]=Level[j-1];
		}

	      else
		{
		  OutliersAws[j]=-1;
		  Level[j]=Level[j-1];
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
	  if (Smoothing[j]!=Smoothing[j-1] && OutliersAws[j-1]==0)
	    {
	      if (j==1 || j==(*l-1))
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
		      if (Smoothing[j]>Smoothing[j-1])
			{                                  
			  OutliersAws[j]=1;
			  Level[j]=Level[j-1];
			}

		      else
			{
			  OutliersAws[j]=-1;
			  Level[j]=Level[j-1];
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

  for (i=1;i<*l;i++)
    {
      if (i==1 || i==(*l-1))
	{
	  if (Region[i]!=Region[i-1])
	    {
	      if (i==1)
		{
		  OutliersAws[0]=1;
		  Region[0]=Region[1];
		}
	      else
		{
		  OutliersAws[*l-1]=1;
		  Region[*l-1]=Region[*l-2];
		}
	    }
	}
      else
	{
	  if (Chromosome[i]!=Chromosome[i-1])
	    {
	      if (Region[i-1]!=Region[i-2])
		{
		  Region[i-1]=Region[i-2];
		  OutliersAws[i-1]=1;
		}

	      if (Region[i+1]!=Region[i])
		{
		  Region[i+1]=Region[i];
		  OutliersAws[i]=1;
		}                
	    }
	  else
	    {
	      if (Region[i]!=Region[i-1] && Region[i+1]!=Region[i] && Region[i+1]==Region[i-1])
		{
		  if (OutliersAws[i-1]==0)
		    {
		      OutliersAws[i]=1;
		      Region[i]=Region[i-1];
		    }
		}
	      else
		{
		  if (Region[i]!=Region[i-1] && OutliersAws[i-1]==0)
		    {
		      Breakpoints[i-1]=1;
		      NextLogRatio[i-1]=LogRatio[i];
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

  for (i=0;i<*l;i++)
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

  *minG=MAXDOUBLE;
  *minAmp=MAXDOUBLE;
  *maxL=-MAXDOUBLE;
  *maxDel=-MAXDOUBLE;

  for (i=0;i<*l;i++)
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

  minG=malloc(1 * sizeof(double));
  maxL=malloc(1 * sizeof(double));
  minAmp=malloc(1 * sizeof(double));
  maxDel=malloc(1 * sizeof(double));

  rangeGainLoss(Smoothing, ZoneGNL, OutliersTot, minG, maxL, minAmp, maxDel, l);


  for (i=0; i<*l;i++)
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
  const int l=*nb_Bkp;
  double UpLeft, LowLeft, UpRight, LowRight, DiffLeft, DiffRight;
  double LRV;
  double SigChr;


  //  for (i in 2:length(profileCGH$BkpInfo[,1]))
  for (i=1;i<l;i++)
    {
      if (BkpInfo_PosOrder[i]==BkpInfo_NextPosOrder[i-1] && BkpInfo_BkpToDel[i-1]==0)
	{
	  SigChr=BkpInfo_Sigma[i];
	  // on regarde d'abord à gauche
	  UpLeft=BkpInfo_Smoothing[i-1] + 3*SigChr;
	  LowLeft=BkpInfo_Smoothing[i-1] - 3*SigChr;

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
		  DiffLeft=fabs(LRV - BkpInfo_Smoothing[i-1]);
		  DiffRight=fabs(LRV - BkpInfo_SmoothingNext[i]);
		  if (DiffRight<DiffLeft)
		    {
		      // On fusionne le Bkp avec la région à droite
		      BkpInfo_BkpToDel[i]=1;
		      BkpInfo_Side[i]=1;
		      BkpInfo_Gap[i-1]=fabs(BkpInfo_Smoothing[i-1]-BkpInfo_SmoothingNext[i]);
		      BkpInfo_Weight[i-1]=1 - kernelpen(BkpInfo_Gap[i-1], *nbsigma*BkpInfo_Sigma[i-1]);
                                
		    }
		  else
		    {
		      // On fusionne le Bkp avec la région à gauche
		      BkpInfo_BkpToDel[i-1]=1;
		      BkpInfo_Side[i-1]=0;
		      BkpInfo_Gap[i]=fabs(BkpInfo_Smoothing[i-1]-BkpInfo_SmoothingNext[i]);
		      BkpInfo_Weight[i]=1 - kernelpen(BkpInfo_Gap[i-1], *nbsigma*BkpInfo_Sigma[i-1]);
                                

		    }
		}
	      else
		{
		  if (((LRV > LowLeft) && (LRV < UpLeft)))
		    {
		      // On fusionne le Bkp avec la région à gauche
		      BkpInfo_BkpToDel[i-1]=1;
		      BkpInfo_Side[i-1]=0;
		      BkpInfo_Gap[i]=fabs(BkpInfo_Smoothing[i-1]-BkpInfo_SmoothingNext[i]);
		      BkpInfo_Weight[i]=1 - kernelpen(BkpInfo_Gap[i-1], *nbsigma*BkpInfo_Sigma[i-1]);
                                
		    }
		  else
		    {
		      // On fusionne le Bkp avec la région à droite
		      BkpInfo_BkpToDel[i]=1;
		      BkpInfo_Side[i]=1;
		      BkpInfo_Gap[i-1]=fabs(BkpInfo_Smoothing[i-1]-BkpInfo_SmoothingNext[i]);
		      BkpInfo_Weight[i-1]=1 - kernelpen(BkpInfo_Gap[i-1], *nbsigma*BkpInfo_Sigma[i-1]);                            
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

void loopTestBkpToMove(double *LogRatio,
		       double *NextLogRatio,
		       double *Smoothing,
		       double *SmoothingNext,
		       int *MoveBkp,
		       int *NbBkp)
{
  int i;
  int l=*NbBkp;

  for (i=0;i<l;i++)
    {
      MoveBkp[i]=testSingle(LogRatio[i],NextLogRatio[i],Smoothing[i],SmoothingNext[i]);
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
