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
	if (it_s_r->second.index.size()>=*msize)
	  {

	    median=median_vector_double(it_s_r->second.LogRatio);
	    mad=mad_vector_double(it_s_r->second.LogRatio);
	    seuil=mad * *alpha;

	    it_LogRatio=it_s_r->second.LogRatio.begin();

	    it_index=it_s_r->second.index.begin();
	    it_index_end=it_s_r->second.index.end();
	    while(it_index != it_index_end)
	      {
		if (*it_LogRatio>median+seuil)
		  {
		    OutliersMad[*it_index]=1;
		    OutliersTot[*it_index]=1;
		  }
		else
		  {
		    if (*it_LogRatio<median-seuil)
		      {
			OutliersMad[*it_index]=-1;
			OutliersTot[*it_index]=-1;
		      }
		  }

		if (OutliersMad[*it_index]==0 & OutliersAws[*it_index]!=0)
		  {
		    OutliersAws[*it_index]=0;
		  }

		if (OutliersMad[*it_index]!=0 & OutliersAws[*it_index]!=0)
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
    int LabelRegionRemoved;
    int LabelRegionRemovedNext;
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
	    agg_aux.Var=var_vector_double(it_LogRatio->second,1);
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





  void putLevel(const double *Smoothing,
		int *Level,
		int *nblevel,
		const int *l)
  {
    int i;
    const int nb=*l;
    map<double, vector<int> > indexLevel;
    map<double, vector<int> >::iterator it_ind;
    map<double, vector<int> >::iterator it_ind_end;
    vector<int>::iterator it_vec;
    vector<int>::const_iterator it_vec_end;

    for (i=0; i<nb; i++)
      {
	indexLevel[Smoothing[i]].push_back(i);
      }

    it_ind=indexLevel.begin();
    it_ind_end=indexLevel.end();

    while(it_ind != it_ind_end)
      {
	*nblevel+=1;
	it_vec=it_ind->second.begin();
	it_vec_end=it_ind->second.end();

	while(it_vec != it_vec_end)
	  {
	    Level[*it_vec]=*nblevel;
	    it_vec++;
	  }

	it_ind++;
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
    int *ZoneGNL=value_dest;
    const double forceGL1=*forceGL1Value;
    const double forceGL2=*forceGL2Value;
    const double NormalRef=*NormalRefValue;
    const double amplicon=*ampliconValue;
    const double deletion=*deletionValue;
    double Smoothing_moins_NormalRef;

    map<int, int > agg_data;

    // construction de la map pour les données aggrégées
    for(i=0;i<*length_src;i++)
      {
	agg_data[index_src[i]]=value_src[i];
      }

    for (i=0;i<*length_dest;i++)
      {
	value_dest[i]=agg_data[index_dest[i]];

	if(NormalRef!=0)
	  {
	    Smoothing_moins_NormalRef=Smoothing[i]-NormalRef;
	  }
	else
	  {
	    Smoothing_moins_NormalRef=Smoothing[i];
	  }

	// Gain et Amplicon
	if(Smoothing_moins_NormalRef>=forceGL2)
	  {
	    if(Smoothing_moins_NormalRef>=amplicon)
	      {
		ZoneGNL[i]=2;
	      }
	    else
	      {
		ZoneGNL[i]=1;
	      }
	  }
	else
	  {
	    if(Smoothing_moins_NormalRef<=forceGL1)
	      {
		if(Smoothing_moins_NormalRef<=deletion)
		  {
		    ZoneGNL[i]=-10;
		  }
		else
		  {
		    ZoneGNL[i]=-1;
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
  double median_vector_double(vector<double> vec)
  {

    sort(vec.begin(),vec.end());
    if(vec.size()%2!=0)
      {
	return(*(vec.begin() + ((vec.size()-1))*0.5));
      }

    else
      {
	return((*(vec.begin() + ((vec.size()))*0.5-1) + *(vec.begin() + ((vec.size()))*0.5))*0.5);
      }

  }

  double mean_vector_double(vector<double> vec)
  {
    return(accumulate(vec.begin(),vec.end(),0.0)/vec.size());

  }

  double mad_vector_double(vector<double> vec)
  {
    double median=0;
    int i;
    const double constant = 1.4826;
    const int vsize=vec.size();
    vector<double> vaux(vsize);
    median=median_vector_double(vec);


    for (i=0;i<vsize;i++)
      {
	vaux[i]=fabs(vec[i]-median);
      }

    return(constant*median_vector_double(vaux));


  }

  double var_vector_double(vector<double> vec, int unbiased)
  {
    double var=0;
    double aux;
    double mean;
    int i;
    const int vsize=vec.size();
    vector<double> vaux(vsize);
    mean=mean_vector_double(vec);

    if (vsize==1)
      return(0);

    for (i=0;i<vsize;i++)
      {
	aux=vec[i]-mean;
	aux*=aux;
	var+=aux;
      }

    if (unbiased==0)
      return(var/(vsize-1));
    else
      return(var/vsize);


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
    if (x>=d)
      return(0);
    else
      {
	k=x/d;
	k=k*k*k;
	k=1-k;
	k=k*k*k;
	return(k);
      }
  }




  double computeSumKernelPen(vector<struct agg> agg_region, double sigma, double d)
  {
    vector<struct agg>::const_iterator it_agg_b = agg_region.begin();
    vector<struct agg>::const_iterator it_agg_n = agg_region.begin();
    vector<struct agg>::const_iterator it_agg_e = agg_region.end();
    double diff=0;
    double sum=0;
    double const inv_sigma=1/sigma;
    ++it_agg_n;
    while(it_agg_n != it_agg_e)
      {
	diff=(*it_agg_n).Mean - (*it_agg_b).Mean;
	diff*=inv_sigma;
	diff=fabs(diff);
	sum+=kernelpen(diff,d);  
	it_agg_b++;
	++it_agg_n;
      }

    return(sum);
  

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
