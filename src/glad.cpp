/*****************************************************************************/
/*                                                                           */
/* Copyright (C) 2004 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2005                            */
/* Contact: bioinfo-staff@curie.fr                                           */
/* It is strictly forbidden to transfer, use or re-use this code             */
/* or part of it without explicit written authorization from Institut Curie. */
/*****************************************************************************/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 
#include <map>
#include <numeric>
#include <functional>
#include <iterator>
#include <math.h>

#ifdef IS_MAC_OS
#include <limits.h>
#else
#include <values.h>
#endif


#ifndef MAXDOUBLE
#include <float.h>
#define MAXDOUBLE DBL_MAX
#endif

#include "glad.h"
#include "glad-function.h"
#include "mva.h"
#include "cutree.h"

#define pi M_PI

using namespace std;

struct paire_double
{
  double value;
  int index_add;
};

struct split_region
{
  vector<double> LogRatio;
  vector<int> index;
};



struct agg
{
  double Mean;
  double Var;
  double VarLike;
  int Card; 
  int LabelRegion;

  agg() {
    Mean = Var = VarLike = 0;
    Card = 0;
  }

  agg(double Mean, double Var, double VarLike, int Card, int LabelRegion) :
    Mean(Mean), Var(Var), VarLike(VarLike), Card(Card), LabelRegion(LabelRegion)
  {
  }

#ifndef USE_COMPARE
  bool operator<(const agg &b) const {
    return LabelRegion < b.LabelRegion;
  }
#endif
};



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
    const int vsize=*l;
    double median;
    double mad;
    double seuil;
    map<int, split_region> s_region;
    map<int, split_region>::iterator it_s_r;
    map<int, split_region>::iterator it_s_r_end;
    vector<int>::iterator it_index;
    vector<int>::iterator it_index_end;
    vector<double>::iterator it_LogRatio;

    for (i=0; i<vsize;i++)
      {
	s_region[Region[i]].LogRatio.push_back(LogRatio[i]);
	s_region[Region[i]].index.push_back(i);
	OutliersTot[i]=0;
	OutliersMad[i]=0;
      }
    it_s_r=s_region.begin();
    it_s_r_end=s_region.end();

    while(it_s_r != it_s_r_end)
      {
	if ((int)(it_s_r->second.index.size()) >= *msize)
	  {

	    median=quantile_vector_double(it_s_r->second.LogRatio, 0.5);
	    mad=mad_vector_double(it_s_r->second.LogRatio);
	    seuil=mad * *alpha;

	    it_LogRatio=it_s_r->second.LogRatio.begin();

	    it_index=it_s_r->second.index.begin();
	    it_index_end=it_s_r->second.index.end();
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

		if ((OutliersMad[*it_index] == 0) & (OutliersAws[*it_index] != 0))
		  {
		    OutliersAws[*it_index]=0;
		  }

		if ((OutliersMad[*it_index] != 0) & (OutliersAws[*it_index] != 0))
		  {
		    OutliersAws[*it_index]=0;
		  }

		if (OutliersAws[*it_index]!=0)
		  {
		    OutliersTot[*it_index]=OutliersAws[*it_index];
		  }


		++it_LogRatio;
		++it_index;
	      }
	  }
	++it_s_r;
      } 
  }
}


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




//   void putLevel(const double *Smoothing,
// 		int *Level,
// 		int *nblevel,
// 		const int *l)
//   {
//     int i;
//     const int nb=*l;
//     map<double, vector<int> > indexLevel;
//     map<double, vector<int> >::iterator it_ind;
//     map<double, vector<int> >::iterator it_ind_end;
//     vector<int>::iterator it_vec;
//     vector<int>::const_iterator it_vec_end;

//     for (i=0; i<nb; i++)
//       {
// 	indexLevel[Smoothing[i]].push_back(i);
//       }

//     it_ind=indexLevel.begin();
//     it_ind_end=indexLevel.end();

//     while(it_ind != it_ind_end)
//       {
// 	*nblevel+=1;
// 	it_vec=it_ind->second.begin();
// 	it_vec_end=it_ind->second.end();

// 	while(it_vec != it_vec_end)
// 	  {
// 	    Level[*it_vec]=*nblevel;
// 	    it_vec++;
// 	  }

// 	it_ind++;
//       }
//   }



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

  void my_merge(const int *index_dest,
		double *value_dest,
		const int *index_src,
		const double *value_src,
		const int *length_dest,
		const int *length_src)
  {
    int i;
    map<int, double > agg_data;

    // construction de la map pour les données aggrégées
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

    // construction de la map pour les données aggrégées
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

    // construction de la map pour les données aggrégées
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

    // construction de la map pour les données aggrégées
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


//   void test_quantile(double *x, const int *n, const double *quantile)
//   {
//     int i;
//     vector<double> vec;
//     double qvalue;

//     for(i = 0; i < *n; i++)
//       {
// 	vec.push_back(x[i]);
//       }

//     qvalue = quantile_vector_double(vec, *quantile);
//     printf("quantile: %f\n", qvalue);

//   }

  double quantile_vector_double(vector<double> vec, const double quantile)
  {
    const double index = (vec.size() - 1) * quantile;
    double h;
    double value_left, value_right;
    double lo, hi;

    lo = floor(index);
    hi = ceil(index);
    h = index - (double)lo;

//     printf("index: %f\n", index);
//     printf("lo: %f\n", lo);
//     printf("hi: %f\n", hi);

//    printf("h: %f\n", h);

    nth_element(vec.begin(),
		vec.begin() + (int)(lo), 
		vec.end());

    value_left = vec[(int)(lo)];
    //    printf("left:%f\n", value_left);



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
//     // de la médiane quantile = 0.5
//     // pour pour quantile = 0.25 et 0.75
//     if(vec.size() % 2)
//       {
// 	// on a un nombre impair d'éléments pour quantile = 0.5
// 	size_t indice = (vec.size() - 1) * quantile;

// 	nth_element(vec.begin(),
// 		    vec.begin() + indice, 
// 		    vec.end());

// 	return vec[indice];
//       }

//     else
//       {
// 	// on a un nombre pair d'éléments pour quantile = 0.5
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
    for(i; i < nb; i++)
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

}




