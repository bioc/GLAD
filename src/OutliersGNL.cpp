/*****************************************************************************/
/* Copyright (C) 2004 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2004                            */
/* Contact: glad@curie.fr                                                    */
/*****************************************************************************/

#include <stdio.h>

#ifdef IS_MAC_OS
#include <limits.h>
#else
#include <values.h>
#endif


#ifndef MAXDOUBLE
#include <float.h>
#define MAXDOUBLE DBL_MAX
#endif

extern "C"
{
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

}
