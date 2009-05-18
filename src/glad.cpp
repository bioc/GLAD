/*****************************************************************************/
/*                                                                           */
/* Copyright (C) 2004 Institut Curie                                         */
/* Author(s): Philippe Hup� (Institut Curie) 2005                            */
/* Contact: bioinfo-staff@curie.fr                                           */
/* It is strictly forbidden to transfer, use or re-use this code             */
/* or part of it without explicit written authorization from Institut Curie. */
/*****************************************************************************/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 
#include <map>
#include <numeric>
#include <functional>
#include <iterator>


#ifdef IS_MAC_OS
#include <limits.h>
#else
#include <values.h>
#endif


#ifndef MAXDOUBLE
#include <float.h>
#define MAXDOUBLE DBL_MAX
#endif

#include "glad-struct.h"
#include "glad.h"
#include "glad-function.h"
#include "mva.h"
#include "cutree.h"


using namespace std;



extern "C" 
{
  void detectOutliers(const double *LogRatio,
		      const int *Region,
		      int *OutliersAws,
		      int *OutliersMad,
		      int *OutliersTot,
		      const int *msize,
		      const double *alpha,
		      const int *l)
  {
    int i;    
    const int vsize = *l;
    double median;
    double mad;
    double seuil;
    map<int, split_region> s_region;
    map<int, split_region>::iterator it_s_r;
    map<int, split_region>::iterator it_s_r_end;
    vector<int>::iterator it_index;
    vector<int>::iterator it_index_end;
    vector<double>::iterator it_LogRatio;

    for (i = 0; i < vsize; i++)
      {
	s_region[Region[i]].LogRatio.push_back(LogRatio[i]);
	s_region[Region[i]].index.push_back(i);
	OutliersTot[i] = 0;
	OutliersMad[i] = 0;
      }
    it_s_r = s_region.begin();
    it_s_r_end = s_region.end();

    while(it_s_r != it_s_r_end)
      {
	if ((int)(it_s_r->second.index.size()) >= *msize)
	  {

	    median = quantile_vector_double(it_s_r->second.LogRatio, 0.5);
	    mad = mad_vector_double(it_s_r->second.LogRatio);
	    seuil = mad * *alpha;

	    it_LogRatio = it_s_r->second.LogRatio.begin();

	    it_index = it_s_r->second.index.begin();
	    it_index_end = it_s_r->second.index.end();
	    while(it_index != it_index_end)
	      {
		if (*it_LogRatio > median + seuil)
		  {
		    OutliersMad[*it_index] = 1;
		    OutliersTot[*it_index] = 1;
		  }
		else
		  {
		    if (*it_LogRatio < median - seuil)
		      {
			OutliersMad[*it_index] = -1;
			OutliersTot[*it_index] = -1;
		      }
		  }

		if ((OutliersMad[*it_index] == 0) && 
		    (OutliersAws[*it_index] != 0))
		  {
		    OutliersAws[*it_index] = 0;
		  }

		if ((OutliersMad[*it_index] != 0) && 
		    (OutliersAws[*it_index] != 0))
		  {
		    OutliersAws[*it_index] = 0;
		  }

		if (OutliersAws[*it_index] != 0)
		  {
		    OutliersTot[*it_index] = OutliersAws[*it_index];
		  }

		++it_LogRatio;
		++it_index;
	      }
	  }
	++it_s_r;
      } 
  }


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
    const int nb=*l;
    int nb_loop=0;
    int stop=0;
    int i,j;
    int nbregion;
    int indRegion;
    double sumkernelpen=0;
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



    while(stop!=1)
      {
	nb_loop++;
	detectOutliers(LogRatio,Region,OutliersAws,OutliersMad,OutliersTot,msize,alpha,l);

	map_ind_reg.erase(map_ind_reg.begin(),map_ind_reg.end());
	LogRatio_no_out.erase(LogRatio_no_out.begin(),LogRatio_no_out.end());
	aggTot.erase(aggTot.begin(),aggTot.end());


	for (i=0;i<nb;i++)
	  {
	    map_ind_reg[Region[i]].push_back(i);

	    if(OutliersTot[i]==0)
	      {
		LogRatio_no_out[Region[i]].push_back(LogRatio[i]);
	      }

	  }

	it_LogRatio=LogRatio_no_out.begin();
	it_LogRatio_end=LogRatio_no_out.end();

	while(it_LogRatio != it_LogRatio_end)
	  {
	    agg_aux.Mean=mean_vector_double(it_LogRatio->second);
	    agg_aux.Var=var_vector_double(it_LogRatio->second, 1);
	    agg_aux.Card=it_LogRatio->second.size();
	    agg_aux.LabelRegion=it_LogRatio->first;

	    if (it_LogRatio->second.size()==1)
	      {
		agg_aux.VarLike=1;
	      }
	    else
	      {
		agg_aux.VarLike=agg_aux.Var;
	      }

	    aggTot.push_back(agg_aux);

	    it_LogRatio++;
	  }

	sort(aggTot.begin(), aggTot.end());


	if (aggTot.size()>1)
	  {

	    sumkernelpen=computeSumKernelPen(aggTot,*sigma, *d);
	    likelihoodGLOBAL=computeLike(aggTot,*lambda, sumkernelpen);

	    nbregion=aggTot.size();


	    min_likelihood=MAXDOUBLE;

	    for(i=0;i<nbregion-1;i++)
	      {
		aggTotAux.erase(aggTotAux.begin(),aggTotAux.end());

		indRegion=i+1;

		for (j=0;j<nbregion;j++)
		  {
		    if(j!=indRegion)
		      {
			aggTotAux.push_back(aggTot[j]);
		      }
		  }

                aggTotAux[i].Card=aggTot[i].Card + aggTot[indRegion].Card;
                barycentre=aggTot[i].Card*aggTot[i].Mean + aggTot[indRegion].Card*aggTot[indRegion].Mean;
                barycentre=barycentre/aggTotAux[i].Card;
                aggTotAux[i].Mean=barycentre;
                within=aggTot[i].Card*aggTot[i].Var + aggTot[indRegion].Card*aggTot[indRegion].Var;
                within=within/aggTotAux[i].Card;
		between=aggTot[i].Card*(aggTot[i].Mean-barycentre)*(aggTot[i].Mean-barycentre);
		between+=aggTot[indRegion].Card*(aggTot[indRegion].Mean-barycentre)*(aggTot[indRegion].Mean-barycentre);
                between=between/aggTotAux[i].Card;
                aggTotAux[i].Var= within + between;
                aggTotAux[i].VarLike=aggTotAux[i].Var;
		sumkernelpen=computeSumKernelPen(aggTotAux,*sigma, *d);

                likelihood=computeLike(aggTotAux,*lambda, sumkernelpen);
                
		if (likelihood<min_likelihood)
		  {
		    min_likelihood=likelihood;
		    LabelRegionRemoved=aggTot[i].LabelRegion;
		    LabelRegionRemovedNext=aggTot[indRegion].LabelRegion;
		  }
	      }

	    if (min_likelihood<likelihoodGLOBAL)
	      {
		it_map_vec=map_ind_reg[LabelRegionRemoved].begin();
		it_map_vec_end=map_ind_reg[LabelRegionRemoved].end();
		while(it_map_vec != it_map_vec_end)
		  {
		    if(Breakpoints[*it_map_vec]==1)
		      {
			Breakpoints[*it_map_vec]=-1;
		      }
		    it_map_vec++;
		  }

		it_map_vec=map_ind_reg[LabelRegionRemovedNext].begin();
		it_map_vec_end=map_ind_reg[LabelRegionRemovedNext].end();
		while(it_map_vec != it_map_vec_end)
		  {
		    Region[*it_map_vec]=LabelRegionRemoved;
		    it_map_vec++;
		  }

	      }
	    else
	      {
		stop=1;
	      }

	  }
	else
	  {
	    stop=1;
	  }
      }
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


    // avec la map, les levels seront ordonn�s par ordre croissant de m�diane
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

  void my_merge(const int *index_dest,
		double *value_dest,
		const int *index_src,
		const double *value_src,
		const int *length_dest,
		const int *length_src)
  {
    int i;
    map<int, double > agg_data;

    // construction de la map pour les donn�es aggr�g�es
    for(i=0;i<*length_src;i++)
      {
	agg_data[index_src[i]]=value_src[i];
      }

    for (i=0;i<*length_dest;i++)
      {
	value_dest[i]=agg_data[index_dest[i]];
      }

  }

  void my_merge_medianlevel(const int *index_dest,
			    int *index_dest_add,
			    double *value_dest,
			    const int *index_src,
			    const int *index_src_add,
			    const double *value_src,
			    const int *length_dest,
			    const int *length_src)
  {
    int i;
    map<int, paire_double > agg_data;

    // construction de la map pour les donn�es aggr�g�es
    for(i=0;i<*length_src;i++)
      {
	agg_data[index_src[i]].value=value_src[i];
	agg_data[index_src[i]].index_add=index_src_add[i];
      }

    for (i=0;i<*length_dest;i++)
      {
	value_dest[i]=agg_data[index_dest[i]].value;
	index_dest_add[i]=agg_data[index_dest[i]].index_add;
      }

  }




  void my_merge_int(const int *index_dest,
		    int *value_dest,
		    const int *index_src,
		    const int *value_src,
		    const int *length_dest,
		    const int *length_src)
  {
    int i;
    map<int, int > agg_data;

    // construction de la map pour les donn�es aggr�g�es
    for(i=0;i<*length_src;i++)
      {
	agg_data[index_src[i]]=value_src[i];
      }

    for (i=0;i<*length_dest;i++)
      {
	value_dest[i]=agg_data[index_dest[i]];
      }

  }

  void my_merge_int_forceGL(const int *index_dest,
			    int *value_dest,
			    const int *index_src,
			    const int *value_src,
			    const int *length_dest,
			    const int *length_src,
			    const double *Smoothing,
			    const double *forceGL1Value,
			    const double *forceGL2Value,
			    const double *NormalRefValue,
			    const double *ampliconValue,
			    const double *deletionValue)
  {
    int i;
    int *ZoneGNL = value_dest;
    const double forceGL1 = *forceGL1Value;
    const double forceGL2 = *forceGL2Value;
    const double NormalRef = *NormalRefValue;
    const double amplicon = *ampliconValue;
    const double deletion = *deletionValue;
    double Smoothing_moins_NormalRef;

    map<int, int > agg_data;

    // construction de la map pour les donn�es aggr�g�es
    for(i = 0; i < *length_src; i++)
      {
	agg_data[index_src[i]] = value_src[i];
      }

    for (i = 0; i < *length_dest; i++)
      {
	value_dest[i] = agg_data[index_dest[i]];

	if(NormalRef!=0)
	  {
	    Smoothing_moins_NormalRef = Smoothing[i] - NormalRef;
	  }
	else
	  {
	    Smoothing_moins_NormalRef = Smoothing[i];
	  }

	// Gain et Amplicon
	if(Smoothing_moins_NormalRef >= forceGL2)
	  {
	    if(Smoothing_moins_NormalRef >= amplicon)
	      {
		ZoneGNL[i] = 2;
	      }
	    else
	      {
		ZoneGNL[i] = 1;
	      }
	  }
	else
	  {
	    if(Smoothing_moins_NormalRef <= forceGL1)
	      {
		if(Smoothing_moins_NormalRef <= deletion)
		  {
		    ZoneGNL[i] = -10;
		  }
		else
		  {
		    ZoneGNL[i] = -1;
		  }
	      }
	  }
      }

  }




  /////////////////////////////////////////////////////////////
  //
  //   Fonctions statistiques
  //
  //////////////////////////////////////////////////////////////

  double IQRdiff(vector<double> vec)
  {
    int i;
    double diffvalue;
    double IQRdiff;
    vector<double> diffvec;

    for (i = 1; i < (int)(vec.size()); i++)
      {
	diffvalue = vec[i] - vec[i - 1];
	diffvec.push_back(diffvalue);
      }

    IQRdiff = IQR_vector_double(diffvec) / 1.908;

    return IQRdiff;
  }

  double IQR_vector_double(vector<double> vec)
  {
    double quantile25;
    double quantile75;

    quantile25 = quantile_vector_double(vec, 0.25);
    quantile75 = quantile_vector_double(vec, 0.75);

    return quantile75 - quantile25;

  }



  double quantile_vector_double(vector<double> vec, const double quantile)
  {
    const double index = (vec.size() - 1) * quantile;
    double h;
    double value_left, value_right;
    double lo, hi;

    lo = floor(index);
    hi = ceil(index);
    h = index - (double)lo;

    nth_element(vec.begin(),
		vec.begin() + (int)(lo), 
		vec.end());

    value_left = vec[(int)(lo)];


    if (h == 0)
      {
	return value_left;
      }
    else
      {
	nth_element(vec.begin(),
		    vec.begin() + (int)(hi), 
		    vec.end());

	value_right = vec[(int)(hi)];
	//	printf("right:%f\n", value_right);

	return (1 - h) * value_left + h * value_right;
      }

  }

//   double quantile_vector_double(vector<double> vec, const double quantile)
//   {

//     // cette fonctionne n'est applicable que pour le calcul
//     // de la m�diane quantile = 0.5
//     // pour pour quantile = 0.25 et 0.75
//     if(vec.size() % 2)
//       {
// 	// on a un nombre impair d'�l�ments pour quantile = 0.5
// 	size_t indice = (vec.size() - 1) * quantile;

// 	nth_element(vec.begin(),
// 		    vec.begin() + indice, 
// 		    vec.end());

// 	return vec[indice];
//       }

//     else
//       {
// 	// on a un nombre pair d'�l�ments pour quantile = 0.5
// 	size_t indice = vec.size() * quantile;
// 	double value_right;
// 	double value_left;

// 	nth_element(vec.begin(),
// 		    vec.begin() + indice - 1, 
// 		    vec.end());

// 	value_left = vec[indice - 1];

// 	nth_element(vec.begin(),
// 		    vec.begin() + indice, 
// 		    vec.end());

// 	value_right = vec[indice];

// 	return (value_left * quantile + value_right * (1 - quantile))  ;
//       }

//   }

  double mean_vector_double(vector<double> vec)
  {
    return accumulate(vec.begin(), vec.end(), 0.0) / vec.size();

  }

  double mad_vector_double(vector<double> vec)
  {
    double median=0;
    int i;
    const double constant = 1.4826;
    const int vsize=vec.size();
    vector<double> vaux(vsize);
    median=quantile_vector_double(vec, 0.5);


    for (i=0;i<vsize;i++)
      {
	vaux[i]=fabs(vec[i]-median);
      }

    return(constant*quantile_vector_double(vaux, 0.5));


  }

  double var_vector_double(vector<double> vec, int unbiased)
  {
    double var = 0;
    double aux;
    double mean;
    int i;
    const int vsize=vec.size();
    vector<double> vaux(vsize);

    mean = mean_vector_double(vec);

    if (vsize == 1)
      return 0;

    for (i = 0; i < vsize; i++)
      {
	aux = vec[i] - mean;
	aux *= aux;
	var += aux;
      }

    if (unbiased == 0)
      return var / (vsize - 1);
    else
      return(var/vsize);


  }


  double median_fabs_double(const double *value, const int l)
  {
    int i;
    vector<double> value_vector;

    for (i = 0; i < l; i++)
      {
	value_vector.push_back(fabs(*value));
	value++;
      }

    return quantile_vector_double(value_vector, 0.5);


  }



  double computeLike(vector<struct agg> agg_region, double lambda, double sumkernelpen)
  {
    vector<struct agg>::iterator it_agg = agg_region.begin();
    vector<struct agg>::iterator it_agg_end = agg_region.end();

    double logsigma=0;
    double nbdata=0;
    while(it_agg != it_agg_end)
      {
	logsigma+=log((*it_agg).VarLike)*(*it_agg).Card;
	nbdata+=(*it_agg).Card;
	it_agg++;
      }
    nbdata=log(nbdata);
    return(logsigma+lambda*sumkernelpen*nbdata);
  }


        
  double kernelpen(double x, const double d)
  {
    double k;
    if (x >= d)
      return 0;
    else
      {
	k = x / d;
	k = k * k * k;
	k = 1 - k;
	k = k * k * k;
	return k;
      }
  }




  double computeSumKernelPen(vector<struct agg> agg_region, double sigma, double d)
  {
    vector<struct agg>::const_iterator it_agg_b = agg_region.begin();
    vector<struct agg>::const_iterator it_agg_n = agg_region.begin();
    vector<struct agg>::const_iterator it_agg_e = agg_region.end();
    double diff = 0;
    double sum = 0;
    double const inv_sigma = 1 / sigma;
    ++it_agg_n;
    while(it_agg_n != it_agg_e)
      {
	diff = (*it_agg_n).Mean - (*it_agg_b).Mean;
	diff *= inv_sigma;
	diff = fabs(diff);
	sum += kernelpen(diff,d);  
	it_agg_b++;
	++it_agg_n;
      }

    return sum;
  
  }



  void printagg(vector<struct agg> agg_region)
  {
    vector<struct agg>::const_iterator b=agg_region.begin();
    vector<struct agg>::const_iterator e=agg_region.end();


    while(b !=e)
      {
	cout << "\tMean=";
	cout <<(*b).Mean;
	cout <<"\tVar=";
	cout <<(*b).Var;
	cout <<"\tVarLike=";
	cout <<(*b).VarLike;
	cout <<"\tCard=";
	cout <<(*b).Card;
	cout <<"\tLabelRegion=";
	cout << (*b).LabelRegion;
	cout << " " << endl;
	b++;
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
				  const int *NbChr,   // Nombre de chromosome � analyser
				  const int *sizeChr, // taille de chaque chromosome
				  const int *startChr, // position pour le d�but des valeurs de chaque chromosome
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

    
    // on prend comme r�f�rence les LogRatios qui sont compris entre certaines + ou - deltaN
    compute_NormalRange(Smoothing,
			NormalRef,
			Level,
			NormalRange,
			deltaN,
			l);

  
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
				   const int *NbChr,   // Nombre de chromosome � analyser
				   const int *sizeChr, // taille de chaque chromosome
				   const int *startChr, // position pour le d�but des valeurs de chaque chromosome
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



}




