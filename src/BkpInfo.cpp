/*****************************************************************************/
/* Copyright (C) 2009 Institut Curie                                         */
/* Author(s): Philippe Hupé (Institut Curie) 2009                            */
/* Contact: glad@curie.fr                                                    */
/*****************************************************************************/

#include "glad-utils.h"

extern "C"
{
 

  void make_BkpInfo(const double *BkpInfo_Gap,
		    int *BkpInfo_GNLchange,
		    const double *BkpInfo_Value,
		    double *BkpInfo_Weight,
		    int *BkpInfo_ZoneGNL,
		    const int *BkpInfo_ZoneGNLnext,
		    const int *nb_Bkp,
		    const double *nbsigma)
  {

    const int l=*nb_Bkp;
    int i;

    for (i=0;i<l;i++)
      {
	BkpInfo_Weight[i]=1 - kernelpen(BkpInfo_Gap[i], *nbsigma*BkpInfo_Value[i]);
	if (BkpInfo_ZoneGNL[i]==BkpInfo_ZoneGNLnext[i])
	  {
	    BkpInfo_GNLchange[i]=0;
	  }
	else
	  {
	    BkpInfo_GNLchange[i]=1;
	  }
      }
  }






}




