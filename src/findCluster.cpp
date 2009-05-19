/*****************************************************************************/
/* Copyright (C) 2009 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2009                            */
/* Contact: glad@curie.fr                                                    */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <map>
#include <numeric>
#include <functional>
#include <algorithm>

#include "glad-struct.h"
#include "glad-utils.h"
#include "findCluster.h"
#include "glad-function.h"
#include "glad.h"
#include "mva.h"
#include "cutree.h"


#define pi M_PI


extern "C"
{

  ///////////////////////
  // Fonction findCluster
  ///////////////////////

  void findCluster(const double *LogRatio,
		   const int *Region,
		   const int *OutliersTot,
		   int *zone,
		   int *method,
		   // paramètres pour clusterglad
		   const double *sigma,
		   const double *d,
		   const double *lambda,
		   const int *nmin,
		   const int *nmax,
		   int *nbclasses,
		   const int *l)
  {

    int i;
    int nc = 1;
    int diag = 0;
    int dist_method = 1;
    int NBR;
    int taille;
    int res = 0;
    int *ia, *ib, *order, *treemerge, *classe, *clusterRegion_Region;
    const int nb = *l;
    double *x_Mean;
    double *dist;
    double *crit;
    double *members;

    map<int, vector<double> > map_Region_LogRatio;
    map<int, vector<double> >::iterator it_map_Region_LogRatio;
    map<int, struct agg> map_clusterRegion;

    // récupération des logratios par région
    for(i = 0; i < nb; i++)
      {
	if(OutliersTot[i] == 0)
	  {
	    map_Region_LogRatio[Region[i]].push_back(LogRatio[i]);
	  }
      }


    NBR = (int)map_Region_LogRatio.size();
    if(NBR == 1)
      {
	*nbclasses = 1;
	for (i = 0; i < nb; i++)
	  {
	    zone[i] = 1;
	  }
      }
    else
      {
	taille = NBR * (NBR - 1) / 2;

	// calcul des informations par région
	x_Mean = (double *)malloc(NBR * sizeof(double));
	members = (double *)malloc(NBR * sizeof(double));
	clusterRegion_Region = (int *)malloc(NBR * sizeof(int));
	it_map_Region_LogRatio = map_Region_LogRatio.begin();

	for (i = 0; i < (int)map_Region_LogRatio.size(); i++)
	  {
	    // on pourrait avoir une fonction qui renvoie Mean et Var pour éviter de calculer 2 fois Mean
	    map_clusterRegion[it_map_Region_LogRatio->first].Mean =  mean_vector_double(it_map_Region_LogRatio->second);
	    x_Mean[i] = map_clusterRegion[it_map_Region_LogRatio->first].Mean;
	    map_clusterRegion[it_map_Region_LogRatio->first].Var = var_vector_double(it_map_Region_LogRatio->second, 1);
	    map_clusterRegion[it_map_Region_LogRatio->first].Card = (int)(it_map_Region_LogRatio->second.size());
	    members[i] = (double)map_clusterRegion[it_map_Region_LogRatio->first].Card;
	    map_clusterRegion[it_map_Region_LogRatio->first].LabelRegion = it_map_Region_LogRatio->first;
	    clusterRegion_Region[i] = it_map_Region_LogRatio->first;

	    // cas des régions avec un seul élément
	    if(it_map_Region_LogRatio->second.size() == 1)
	      {
		map_clusterRegion[it_map_Region_LogRatio->first].Var = 0;
		map_clusterRegion[it_map_Region_LogRatio->first].VarLike = 1;
	      }
	    else
	      {
		map_clusterRegion[it_map_Region_LogRatio->first].VarLike = map_clusterRegion[it_map_Region_LogRatio->first].Var;
	      }

	    it_map_Region_LogRatio++;
	  }

	// calcul de la matrice de distance
	dist = (double *)calloc(taille, sizeof(double));

	R_distance(x_Mean,
		   &NBR,
		   &nc, 
		   dist, 
		   &diag , 
		   &dist_method);

	free(x_Mean);

	// clustering hiérarchique
	ia = (int *)calloc(NBR, sizeof(int));
	ib = (int *)calloc(NBR, sizeof(int));
	order = (int *)calloc(NBR, sizeof(int));
	crit = (double *)calloc(NBR, sizeof(double));


	hclust(&NBR,
	       &taille, 
	       method,
	       ia,
	       ib,
	       order,
	       crit,
	       members,
	       dist,
	       &res);


	free(dist);
	free(members);
	free(order);
	free(crit);
 
	treemerge = (int *)malloc((2 * NBR - 2) * sizeof(int));
	memcpy(&treemerge[0], &ia[0], (NBR - 1) * sizeof(int));
	memcpy(&treemerge[NBR - 1], &ib[0], (NBR - 1) * sizeof(int));

	free(ia);
	free(ib);


	// calcul du nombre de classes
	*nbclasses = clusterglad(map_clusterRegion,
				 treemerge,
				 *nmin,
				 *nmax,
				 *sigma,
				 *d,
				 *lambda);

	classe = (int *)malloc(NBR * sizeof(int));

	R_cutree(treemerge,
		 nbclasses,
		 classe,
		 &NBR);

	free(treemerge);

	my_merge_int(Region,
		     zone,
		     clusterRegion_Region,
		     classe,
		     l,
		     &NBR);

	free(classe);
	free(clusterRegion_Region);
      }
  }


  int clusterglad(map<int, struct agg> map_clusterRegion,
		  int *treemerge,
		  const int min,
		  const int max,
		  const double sigma,
		  const double d,
		  const double lambda)
  {
    int i, j;
    int NBR = (int)map_clusterRegion.size();
    int nmin = min;
    int nmax = max;
    int NbTotObs = 0;
    int *classe;
    double logVar;
    double Mean;
    double logLike;
    double min_logLike;
    double sum_kernelpen;
    vector<double> vec_logVar;
    vector<double> vec_Mean;
    vector<double> vec_logLike;
    vector<double> vec_deltaoversigma;

    map<int, struct agg>::iterator it_map_clusterRegion;


  if(nmin > nmax)
    {
      printf("in clusterglad function: nmin greater than nmax\n");
    }

  if (nmin == nmax)
    {
      if(NBR > nmin)
	{
	  return nmin;
	}
      else
	{
	  return NBR;
	}
    }



  if (nmax > NBR)
    {
      nmax = NBR;
    }

  it_map_clusterRegion = map_clusterRegion.begin();
  for (i =0; i < NBR; i++)
    {
      NbTotObs += it_map_clusterRegion->second.Card;
      it_map_clusterRegion++;
    }

  //  printf("NbTotObs: %i\n", NbTotObs);


  classe = (int *) malloc(NBR * sizeof(int));
  for (i = nmin; i <= nmax; i++)
    {


      R_cutree(treemerge,
	       &i,
	       classe,
	       &NBR);

//       printf("classe[%i]\n", i);
//       for (j = 0; j < NBR; j++)
// 	{
// 	  printf("%i\t", classe[j]);
// 	}
//       printf("\n");

      for (j = 1; j <= i; j++)
	{
	  mergeLike(map_clusterRegion,
		    &logVar,
		    &Mean,
		    classe,
		    j);

	  vec_logVar.push_back(logVar);
	  vec_Mean.push_back(Mean);
	}

//       for (j = 0; j < (int)vec_logVar.size(); j++)
// 	{
// 	  printf("logVar: %f - Mean: %f\n", vec_logVar[j], vec_Mean[j]);
// 	}

      sort(vec_Mean.begin(), vec_Mean.end());

      // calcul de deltaoversigma
      for (j = 1; j < (int)vec_Mean.size(); j++)
	{
	  vec_deltaoversigma.push_back(fabs(vec_Mean[j] - vec_Mean[j - 1]) / sigma);
	  //	  printf("deltaoversigma: %f\n", vec_deltaoversigma[j - 1]);
	}

      // calcul de sum_kernelpen
      if((int)vec_Mean.size() == 1)
	{
	  sum_kernelpen = 1;
	}
      else
	{
	  sum_kernelpen = 0;
	}

      for (j = 0; j < (int)vec_deltaoversigma.size(); j++)
	{
	  sum_kernelpen += kernelpen(vec_deltaoversigma[j], d);
	}

      logLike = 0;
      logLike += lambda * sum_kernelpen * log(NbTotObs);
      logLike += accumulate(vec_logVar.begin(), vec_logVar.end(), 0.0);


      vec_logLike.push_back(logLike);

      //      printf("logLike[%i]: %f\n", i, logLike);


      vec_logVar.erase(vec_logVar.begin(), vec_logVar.end());
      vec_Mean.erase(vec_Mean.begin(), vec_Mean.end());
      vec_deltaoversigma.erase(vec_deltaoversigma.begin(), vec_deltaoversigma.end());

    }

  free(classe);


  min_logLike = *min_element(vec_logLike.begin(), vec_logLike.end());

  //  printf("minlogLike: %f\n", min_logLike);
  //  return(nmin + which(logLike == min(logLike))[1] - 1)

  for (i = 0; i < (int)vec_logLike.size(); i++)
    {
      if(min_logLike == vec_logLike[i])
	{
	  return nmin + i;
	}
    }

  return 0;
  }

  void mergeLike(map<int, struct agg> map_clusterRegion,
		 double *logVar,
		 double *Mean,
		 const int *classe,
		 const int which_classe)
  {
    int nbobs = 0;
    double barycentre = 0;
    double between = 0;
    double tmp_between;
    double within = 0;
    double variance;
    int i;
    map<int, struct agg>::iterator it_map_clusterRegion;

    it_map_clusterRegion = map_clusterRegion.begin();
    for (i = 0; i < (int)map_clusterRegion.size(); i++)
      {
	if (classe[i] == which_classe)
	  {
	    nbobs += it_map_clusterRegion->second.Card;
	    barycentre += it_map_clusterRegion->second.Card * it_map_clusterRegion->second.Mean;
	    within += it_map_clusterRegion->second.Card * it_map_clusterRegion->second.Var;
	    it_map_clusterRegion++;
	  }
	else
	  {
	    it_map_clusterRegion++;
	  }
      }

    barycentre /= nbobs;
    within /= nbobs;

    it_map_clusterRegion = map_clusterRegion.begin();
    for (i = 0; i < (int)map_clusterRegion.size(); i++)
      {
	if (classe[i] == which_classe)
	  {
	    tmp_between = it_map_clusterRegion->second.Mean - barycentre;
	    tmp_between = tmp_between * tmp_between;
	    between += it_map_clusterRegion->second.Card * tmp_between;
	    it_map_clusterRegion++;
	  }
	else
	  {
	    it_map_clusterRegion++;
	  }
      }

    between /= nbobs;

    variance = within + between;


    if (nbobs == 0)
      {
	*logVar = 0;
	*Mean = barycentre;
      }
    else
      {
	*logVar = nbobs * (log(variance) + (1 + log(2 * pi)));
	*Mean = barycentre;
      }
  }  

  void filterBkp_moveBkp_Outliers (int *ZoneGNL,
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
    const int nb = *l;
    const int nb_moins_un = *l - 1;
    const int nb_moins_deux= *l - 2;

    for (pos = 1; pos < nb; pos++)
      {
	pos_moins_un = pos - 1;
	if (Chromosome[pos] == Chromosome[pos_moins_un])
	  {
	    pos_plus_un = pos + 1;
	    if (OutliersTot[pos] != 0 && 
		Breakpoints[pos] == 1 && 
		ZoneGNL[pos] == ZoneGNL[pos_plus_un] && 
		ZoneGNL[pos_moins_un] != ZoneGNL[pos_plus_un])
	      {
		*RecomputeSmt = 1;
		Breakpoints[pos] = 0;
		Breakpoints[pos_moins_un] = 1;
		OutliersTot[pos] = 0;
		OutliersAws[pos] = 0;
		Level[pos] = Level[pos_plus_un];

	      }

	    if (pos < nb_moins_un && 
		Breakpoints[pos] == 1 && 
		OutliersTot[pos_plus_un] != 0 && 
		ZoneGNL[pos] == ZoneGNL[pos_plus_un] && 
		ZoneGNL[pos_moins_un] != ZoneGNL[pos_plus_un])
	      {
		*RecomputeSmt = 1;	  
		Breakpoints[pos] = 0;
		Breakpoints[pos_plus_un] = 1;
		OutliersTot[pos_plus_un] = 0;
		OutliersAws[pos_plus_un] = 0;
		Level[pos_plus_un] = Level[pos];

		if (pos < nb_moins_deux)
		  {
		    pos_plus_deux = pos + 2;
		    if (Chromosome[pos_plus_un] == Chromosome[pos_plus_deux])
		      {
			if (Level[pos_plus_un] == Level[pos_plus_deux])
			  {
			    Breakpoints[pos_plus_un] = 0;
			  }
		      }		    		        
		  }
	      }
	  }
      }
  }
}

