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
#include <values.h>

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
                  int *l)
{
  int pos;
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
  for (pos=1;pos<(*l-1);pos++)
    {

      if (Level[pos-1]==Level[pos+1] & Level[pos-1]!=Level[pos])
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
      
      if (Level[pos-1]==Level[pos+1] & Level[pos-1]!=Level[pos])
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
	  if (OutliersTot[pos]!=0 & Breakpoints[pos]==1 & ZoneGNL[pos]==ZoneGNL[pos+1] & ZoneGNL[pos-1]!=ZoneGNL[pos+1])
	    {
	      *RecomputeSmt=1;
	      Breakpoints[pos]=0;
	      Breakpoints[pos-1]=1;
	      OutliersTot[pos]=0;
	      OutliersAws[pos]=0;
	      Level[pos]=Level[pos+1];

	    }

	  if (pos < (*l-1) & Breakpoints[pos]==1 & OutliersTot[pos+1]!=0 & ZoneGNL[pos]==ZoneGNL[pos+1] & ZoneGNL[pos-1]!=ZoneGNL[pos+1])
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

      if (Smoothing[j]!=Smoothing[j-1] & Smoothing[j+1]!=Smoothing[j] & Smoothing[j+1]==Smoothing[j-1] & j<(*l-1))
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
	  if (Smoothing[j]!=Smoothing[j-1] & OutliersAws[j-1]==0)
	    {
	      if (j==1 | j==(*l-1))
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
      if (i==1 | i==(*l-1))
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
	      if (Region[i]!=Region[i-1] & Region[i+1]!=Region[i] & Region[i+1]==Region[i-1])
		{
		  if (OutliersAws[i-1]==0)
		    {
		      OutliersAws[i]=1;
		      Region[i]=Region[i-1];
		    }
		}
	      else
		{
		  if (Region[i]!=Region[i-1] & OutliersAws[i-1]==0)
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

/* 	  if (ZoneGNL[i]==1) */
/* 	    { */
/* 	      if (Smoothing[i]<*minG) */
/* 		{ */
/* 		  *minG=Smoothing[i]; */
/* 		} */
/* 	    } */

/* 	  if(ZoneGNL[i]==-1) */
/* 	    { */
/* 	      if (Smoothing[i]>*maxL) */
/* 		{ */
/* 		  *maxL=Smoothing[i]; */
/* 		} */
/* 	    } */

/* 	  if (ZoneGNL[i]==2) */
/* 	    { */
/* 	      if(Smoothing[i]<*minAmp) */
/* 		{ */
/* 		  *minAmp=Smoothing[i]; */
/* 		} */
/* 	    } */

/* 	  if (ZoneGNL[i]==-10) */
/* 	    { */
/* 	      if (Smoothing[i]>*maxDel) */
/* 		{ */
/* 		  *maxDel=Smoothing[i]; */
/* 		} */
/* 	    } */
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

      if (Smoothing[i]>=*minG | Smoothing[i]>=*minAmp)
	{
	  ZoneGNL[i]=1;
	 
	  if(Smoothing[i]>=*minAmp)
	    {
	      ZoneGNL[i]=2;
	    }

	}

      else
	{
	  if (Smoothing[i]<=*maxL | Smoothing[i]<=*maxDel)
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
