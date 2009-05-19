/*****************************************************************************/
/* Copyright (C) 2009 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2009                            */
/* Contact: glad@curie.fr                                                    */
/*****************************************************************************/

#include <map>
#include <vector>

#include "glad-struct.h"
#include "loopRemove.h"

extern "C"
{
  void loopRemove(const double *LogRatio,
		  int *Region,
		  int *OutliersAws,
		  int *OutliersMad,
		  int *OutliersTot,
		  int *Breakpoints,
		  const int *msize,
		  const double *alpha,
		  const double *lambda,
		  const double *d,
		  const double *sigma,
		  const int *l)
  {
    const int nb = *l;
    int nb_loop = 0;
    int stop = 0;
    int i, j;
    int nbregion;
    int indRegion;
    double sumkernelpen = 0;
    double likelihoodGLOBAL, likelihood;
    double barycentre, within, between;
    double min_likelihood;
    int LabelRegionRemoved = 0;
    int LabelRegionRemovedNext = 0;
    map<int, vector<int> > map_ind_reg;
    vector<int>::iterator it_map_vec;
    vector<int>::const_iterator it_map_vec_end;
    vector<struct agg> aggTot;
    vector<struct agg> aggTotAux;
    struct agg agg_aux;
    map<int, vector<double> > LogRatio_no_out;
    map<int, vector<double> >::const_iterator it_LogRatio;
    map<int, vector<double> >::const_iterator it_LogRatio_end;



    while(stop != 1)
      {
	nb_loop++;
	detectOutliers(LogRatio,Region,OutliersAws,OutliersMad,OutliersTot,msize,alpha,l);

	map_ind_reg.erase(map_ind_reg.begin(),map_ind_reg.end());
	LogRatio_no_out.erase(LogRatio_no_out.begin(),LogRatio_no_out.end());
	aggTot.erase(aggTot.begin(),aggTot.end());


	for (i = 0; i < nb; i++)
	  {
	    map_ind_reg[Region[i]].push_back(i);

	    if(OutliersTot[i] == 0)
	      {
		LogRatio_no_out[Region[i]].push_back(LogRatio[i]);
	      }

	  }

	it_LogRatio = LogRatio_no_out.begin();
	it_LogRatio_end = LogRatio_no_out.end();

	while(it_LogRatio != it_LogRatio_end)
	  {
	    agg_aux.Mean = mean_vector_double(it_LogRatio->second);
	    agg_aux.Var = var_vector_double(it_LogRatio->second, 1);
	    agg_aux.Card = it_LogRatio->second.size();
	    agg_aux.LabelRegion = it_LogRatio->first;

	    if (it_LogRatio->second.size() == 1)
	      {
		agg_aux.VarLike = 1;
	      }
	    else
	      {
		agg_aux.VarLike = agg_aux.Var;
	      }

	    aggTot.push_back(agg_aux);

	    it_LogRatio++;
	  }

	sort(aggTot.begin(), aggTot.end());


	if (aggTot.size()>1)
	  {

	    sumkernelpen = computeSumKernelPen(aggTot,*sigma, *d);
	    likelihoodGLOBAL = computeLike(aggTot,*lambda, sumkernelpen);

	    nbregion = aggTot.size();


	    min_likelihood = MAXDOUBLE;

	    for(i = 0; i< (nbregion - 1); i++)
	      {
		aggTotAux.erase(aggTotAux.begin(),aggTotAux.end());

		indRegion = i + 1;

		for (j = 0; j < nbregion; j++)
		  {
		    if(j != indRegion)
		      {
			aggTotAux.push_back(aggTot[j]);
		      }
		  }

                aggTotAux[i].Card = aggTot[i].Card + aggTot[indRegion].Card;
                barycentre = aggTot[i].Card*aggTot[i].Mean + aggTot[indRegion].Card*aggTot[indRegion].Mean;
                barycentre = barycentre / aggTotAux[i].Card;
                aggTotAux[i].Mean = barycentre;
                within = aggTot[i].Card * aggTot[i].Var + aggTot[indRegion].Card * aggTot[indRegion].Var;
                within = within / aggTotAux[i].Card;
		between = aggTot[i].Card * (aggTot[i].Mean-barycentre) * (aggTot[i].Mean-barycentre);
		between += aggTot[indRegion].Card * (aggTot[indRegion].Mean-barycentre) * (aggTot[indRegion].Mean-barycentre);
                between = between / aggTotAux[i].Card;
                aggTotAux[i].Var = within + between;
                aggTotAux[i].VarLike=aggTotAux[i].Var;
		sumkernelpen=computeSumKernelPen(aggTotAux,*sigma, *d);

                likelihood = computeLike(aggTotAux,*lambda, sumkernelpen);
                
		if (likelihood < min_likelihood)
		  {
		    min_likelihood = likelihood;
		    LabelRegionRemoved = aggTot[i].LabelRegion;
		    LabelRegionRemovedNext = aggTot[indRegion].LabelRegion;
		  }
	      }

	    if (min_likelihood < likelihoodGLOBAL)
	      {
		it_map_vec = map_ind_reg[LabelRegionRemoved].begin();
		it_map_vec_end = map_ind_reg[LabelRegionRemoved].end();
		while(it_map_vec != it_map_vec_end)
		  {
		    if(Breakpoints[*it_map_vec] == 1)
		      {
			Breakpoints[*it_map_vec] = -1;
		      }
		    it_map_vec++;
		  }

		it_map_vec = map_ind_reg[LabelRegionRemovedNext].begin();
		it_map_vec_end = map_ind_reg[LabelRegionRemovedNext].end();
		while(it_map_vec != it_map_vec_end)
		  {
		    Region[*it_map_vec] = LabelRegionRemoved;
		    it_map_vec++;
		  }

	      }
	    else
	      {
		stop = 1;
	      }

	  }
	else
	  {
	    stop = 1;
	  }
      }
  }


  ///////////////////////
  // Fonction removeLevel
  ///////////////////////

  void loop_chromosome_removeLevel(const double *LogRatio,
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
				   const int *BkpDetected)


  {

    int i;
    int start, size;



    for(i = 0; i < *NbChr; i++)
      {
	start = startChr[i];
	size = sizeChr[i];


	if(!BkpDetected[i])
	  {
	    detectOutliers(&LogRatio[start],
			   &Level[start],
			   &OutliersAws[start],
			   &OutliersMad[start],
			   &OutliersTot[start],
			   msize,
			   alpha,
			   &size);
	  }
	else
	  {
	    // Optimisation des Breakpoints
	    loopRemove(&LogRatio[start],
		       &Level[start],
		       &OutliersAws[start],
		       &OutliersMad[start],
		       &OutliersTot[start],
		       &Breakpoints[start],
		       msize,
		       alpha,
		       lambda,
		       d,
		       &sigma[i],
		       &size);

	    updateBkpRL(&Level[start],
			&OutliersAws[start],
			&Breakpoints[start],
			&PosOrder[start],
			&NextLogRatio[start],
			&LogRatio[start],
			&size);

	    detectOutliers(&LogRatio[start],
			   &Level[start],
			   &OutliersAws[start],
			   &OutliersMad[start],
			   &OutliersTot[start],
			   msize,
			   alpha,
			   &size);
	  }
      }
  }

  /////////////////////////////////////////////////
  // Fonction Optimisation du nombre de Breakpoints
  ////////////////////////////////////////////////


  void OptmisationBreakpointsStep(double *Smoothing,
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
				  const int *l) // nombre total de sondes
  {
    loop_chromosome_removeLevel(LogRatio,
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
				NbChr,   
				sizeChr,
				startChr,
				BkpDetected);

    // calcul de la smoothing line
    compute_median_smoothing(LogRatio,
			     Level,
			     Smoothing,
			     l);

    
    // on prend comme référence les LogRatios qui sont compris entre certaines + ou - deltaN
    compute_NormalRange(Smoothing,
			NormalRef,
			Level,
			NormalRange,
			deltaN,
			l);

  
  }


}
