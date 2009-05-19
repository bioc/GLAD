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
#include <math.h>

#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 
#include <map>
#include <numeric>
#include <functional>
#include <iterator>



#include "glad-struct.h"
#include "glad-utils.h"
#include "glad-function.h"


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


}
