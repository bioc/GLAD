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


/*************************************/
/* fonctions utilisées dans filterBkp */
/*************************************/

void updateLevel (int *Chromosome,
                  int *Breakpoints,
                  int *Level,
                  int *NextPosOrder,
                  int *PosOrder,
		  double *NextLogRatio,
		  double *LogRatio,
		  int *BeforePosOrder,
                  int *l)
{
  /*printf("External Call: updateLevel\n");*/
  int pos;
  /*int comp;*/
  for (pos=1;pos<*l;pos++)
    {
      /*comp=strcmp(Chromosome[pos],Chromosome[pos-1]);*/
      if (Chromosome[pos]==Chromosome[pos-1])
	{
	  if (Breakpoints[pos-1]!=1)
	    {
	      Level[pos]=Level[pos-1];
	    }
	  if (Breakpoints[pos-1]==1)
	    {
	      NextPosOrder[pos-1]=PosOrder[pos];
	      NextLogRatio[pos-1]=LogRatio[pos];
	      if (pos>1)
		{
		  BeforePosOrder[pos-1]=PosOrder[pos-2];
		}
	    }
	}
    }

}



void updateOutliers (int *OutliersAws,
                     int *Level,
                     int *Breakpoints,
                     int *l)
{
  /*printf("External Call: updateOutliers\n");*/
  
  int pos;
  for (pos=1;pos<(*l-1);pos++)
    {

      if (Level[pos-1]==Level[pos+1] & Level[pos-1]!=Level[pos])
	{
	  Level[pos]=Level[pos-1];
	  Breakpoints[pos-1]=0;
	  Breakpoints[pos]=0;
	  OutliersAws[pos]=1;
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
  /*printf("External Call: moveBkp\n");*/

  int pos;
  /*int comp1;*/
  /*int comp2;*/

  for (pos=1;pos<*l;pos++)
    {
      /*comp1=strcmp(Chromosome[pos],Chromosome[pos-1]);*/
      if (Chromosome[pos]==Chromosome[pos-1])
	{
	  if (OutliersTot[pos]!=0 & Breakpoints[pos]==1 & ZoneGNL[pos]==ZoneGNL[pos+1] & ZoneGNL[pos-1]!=ZoneGNL[pos+1])
	    {
	      printf("gauche moveBkp: %i\n",pos);
	      *RecomputeSmt=1;
	      Breakpoints[pos]=0;
	      Breakpoints[pos-1]=1;
	      OutliersTot[pos]=0;
	      OutliersAws[pos]=0;
	      Level[pos]=Level[pos+1];

	    }

	  if (pos < (*l-1) & Breakpoints[pos]==1 & OutliersTot[pos+1]!=0 & ZoneGNL[pos]==ZoneGNL[pos+1] & ZoneGNL[pos-1]!=ZoneGNL[pos+1])
	    {
	      printf("droite moveBkp: %i\n",pos);
	      *RecomputeSmt=1;	  
	      Breakpoints[pos]=0;
	      Breakpoints[pos+1]=1;
	      OutliersTot[pos+1]=0;
	      OutliersAws[pos+1]=0;
	      Level[pos+1]=Level[pos];

	      if (pos < (*l-2))
		{
		  /*comp2=strcmp(Chromosome[pos+1],Chromosome[pos+2]);*/

		  if (Chromosome[pos+1]==Chromosome[pos+2])
		    {
		      if (Level[pos+1]==Level[pos+2])
			{
			  /*printf("droite: il faut supprimer le Bkp\n");*/
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
  /*printf("External Call: awsBkp\n");*/

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
                  int *NextPosOrder,
		  double *NextLogRatio,
		  double *LogRatio,
		  int *BeforePosOrder,
                  int *l)
{
  /*printf("External Call: updateBkpRL\n");*/

  int i;
  /*int comp;*/

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
	  /*comp=strcmp(Chromosome[i],Chromosome[i-1]);*/
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
		      NextPosOrder[i-1]=PosOrder[i];
		      NextLogRatio[i-1]=LogRatio[i];
		      BeforePosOrder[i-1]=PosOrder[i-2];
		    }
		}                                               
	    }
	}
    }
}
