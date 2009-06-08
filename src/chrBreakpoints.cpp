/*****************************************************************************/
/* Copyright (C) 2009 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2009                            */
/* Contact: glad@curie.fr                                                    */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <map>
#include <vector>
#include <algorithm>

using namespace std;

#include "glad-limits.h"
#include "glad-struct.h"
#include "chrBreakpoints.h"
#include "glad-utils.h"
#include "HaarSeg.h"



extern "C"
{

  /*******************************************/
  /* fonctions utilisées dans chrBreakpoints */
  /*******************************************/
  void chrBreakpoints_haarseg(const double *LogRatio,
			      const int *Chromosome,
			      double *Smoothing,
			      int *Level,
			      int *OutliersAws,
			      int *regionChr,
			      int *Breakpoints,
			      int *sizeChr, // taille de chaque chromosome
			      int *startChr, // position pour le début des valeurs de chaque chromosome
			      int *IQRChr, // numéro du chromosome pour le calcul de l'IQR
			      double *IQRValue, // valeur de l'IQR
			      int *BkpDetected,
			      // paramètres pour Haarseg
			      const double *breaksFdrQ,
			      const int *haarStartLevel,
			      const int *haarEndLevel,
			      const int *NbChr, // nombre de chromosomes
			      const int *l, // nombre de probes
			      double *weights) // poids de chacune des probes
  {
    int i;
    int start, size;
    int nblevel = 0;
    int nbregion = 0;
    int stepHalfSize;
    int *peakLoc;
    const int nb = *l;
    double *convResult;
    double *weights_aux = NULL;

    map<int, vector<double> > LogRatio_byChr;
    map<int, vector<double> >::iterator it_LogRatio_byChr;


    // récupération des informations sur la taille des chromosomes
    // et calcul de la variance par chromosome
    for (i = 0; i < nb ; i++)
      {
	LogRatio_byChr[Chromosome[i]].push_back(LogRatio[i]);
      }

    startChr[0] = 0;
    it_LogRatio_byChr = LogRatio_byChr.begin();
    for (i = 0; i < *NbChr; i++)
      {
	sizeChr[i] = it_LogRatio_byChr->second.size();
	IQRChr[i] = it_LogRatio_byChr->first;
	IQRValue[i] = IQRdiff(it_LogRatio_byChr->second);

	if(i > 0)
	  {
	    startChr[i] = startChr[i - 1] + sizeChr[i - 1];
	    if (IQRChr[i] < IQRChr[i - 1])
	      {
		printf("WARNINGS: Chromosome are not correctly ordered\n");
		printf("i:%i - i+1:%i\n", IQRChr[i], IQRChr[i - 1]);
	      }
	  }
	it_LogRatio_byChr++;
      }

    // segmentation chromosome par chromosome
    for (i = 0; i < *NbChr; i++)
      {
	start = startChr[i];
	size = sizeChr[i];
	stepHalfSize = 1;

	convResult = (double *)calloc(size, sizeof(double));
	peakLoc = (int *)calloc(size, sizeof(int));


	printf("Chromosome %i\n", i + 1);
	if(weights != NULL)
	  {
	    weights_aux = &weights[start];
	  }

	HaarSegGLAD(&LogRatio[start],
		    &size,
		    &stepHalfSize,
		    convResult,
		    peakLoc,
                    breaksFdrQ,
                    haarStartLevel,
                    haarEndLevel,
		    &Smoothing[start],
		    weights_aux);


	free(convResult);
	free(peakLoc);

	nbregion += 1;

	putLevel_awsBkp(&Smoothing[start],
                        &LogRatio[start],                              
			&Level[start],  
			&nblevel,       
			&size,
			&OutliersAws[start], 
			&nbregion,           
			&regionChr[start],   
			&Breakpoints[start], 
			&BkpDetected[i]);     

      }

  }


  void putLevel_awsBkp(//variables pour putLevel
		       double *Smoothing,
		       const double *LogRatio,
		       int *Level,
		       int *nblevel,
		       const int *l,
		       // variables pour awsBkp
		       int *OutliersAws,
		       int *nbregion,
		       int *regionChr,
		       int *Breakpoints,
		       int *bkp_detected)
  {

    putLevel(Smoothing,
	     LogRatio,
	     Level,
	     nblevel,
	     l);

    awsBkp (Smoothing,
	    OutliersAws,
	    Level,
	    nbregion,
	    regionChr,
	    Breakpoints,
	    bkp_detected,
	    l);


  }

  void awsBkp(const double *Smoothing,
	      int *OutliersAws,
	      int *Level,
	      int *nbregion,
	      int *regionChr,
	      int *Breakpoints,
	      int *bkp_detected,
	      const int *l)
  {

    int j;
    int j_moins_un;
    int j_plus_un;
    int last_bkp_pos = -1;

    int *rupture;

    const int nb = *l;
    const int nb_moins_un = *l - 1;

    rupture = (int *)calloc(nb, sizeof(int));

    // initialisation de la première valeur de regionChr à nbregion
    regionChr[0] = *nbregion;

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

    memcpy(Breakpoints, &rupture[1], (nb - 1) * sizeof(int)); // <-> subsetdata$Breakpoints <- c(awsBkp$rupture[2:l],0)

    free(rupture);
  }

  void putLevel(double *Smoothing,
		const double *LogRatio,
		int *Level,
		int *nblevel,
		const int *l)
  {
    int i, j;
    const int nb = *l;
    double MedianValue;
    double SmoothingValue;

    map<double, vector<double> > LogRatioLevel;
    map<double, vector<double> >::iterator it_LogRatioLevel;
    map<double, double> MedianLevel;
    map<double, double>::iterator it_MedianLevel;
    map<double, vector<int> > indexLevel;
    vector<int>::iterator it_vec;

    for (i = 0; i < nb; i++)
      {
	indexLevel[Smoothing[i]].push_back(i);
	LogRatioLevel[Smoothing[i]].push_back(LogRatio[i]);
      }


    // avec la map, les levels seront ordonnés par ordre croissant de médiane
    it_LogRatioLevel = LogRatioLevel.begin();
    for (i = 0; i < (int)LogRatioLevel.size(); i++)
      {
	MedianLevel[it_LogRatioLevel->first] = quantile_vector_double(it_LogRatioLevel->second, 0.5);
	it_LogRatioLevel++;
      }


    it_MedianLevel = MedianLevel.begin();
    for(i = 0; i < (int)MedianLevel.size(); i++)
      {
	*nblevel += 1;

	SmoothingValue = it_MedianLevel->first;
	MedianValue = it_MedianLevel->second;

	it_vec = indexLevel[SmoothingValue].begin();
	for(j = 0; j < (int)(indexLevel[SmoothingValue].size()); j++)
	  {
	    Level[*it_vec] = *nblevel;
	    Smoothing[*it_vec] = MedianValue;
	    it_vec++;
	  }

	it_MedianLevel++;
      }

  }


}
