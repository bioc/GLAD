#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "laws.h"

#define max(a,b) ((a)>=(b)?(a):(b))
#define min(a,b) ((a)<=(b)?(a):(b))

/*Le code est fait pour une fonction gaussienne avec l'option symmetric=TRUE*/



void iawsuni (double *y,
              int *n,
	      double *hinit,
	      double *bi,
	      double *ai,
	      double *kern)
{

  int i,j;
  int ja,je,ih,iz;
  double z,wj,az,swj,swjy;
  double diff;
  double inv_hinit=10/(*hinit);


  ih=*hinit;


  for (i=1;i<=*n;i++)
    {
      ja=max(1,i-ih);
      je=min(*n,i+ih);
      swj=0;
      swjy=0;
      for (j=ja;j<=je;j++)
	{
	  diff=i-j;	  
	  z=diff*inv_hinit;
	  /*z=z*z*100;*/
	  /* Comme inv_hinit=*hinit*10 on a: */
	  z=z*z;
	  /*printf("z vaut %.9f\n",z);*/
	  if (z<100)
	    {	      	   
	      iz=z;
	      az=z-iz;
	      wj=(kern[iz+1] - kern[iz])*az + kern[iz];
	      swj=swj+wj;
	      swjy=swjy+wj*y[j-1];	   
	    }	
	}

      ai[i-1]=swjy;
      bi[i-1]=swj;
    }

}   




void lawsuni(double *y,
	     const int n,
	     double inv_hakt,
	     int ih,
	     const double inv_lamakt,
	     double *theta,
	     double *bi,
	     double *ai,
	     double *kernl,
	     double *kerns)
{
  int i,j;
  int ja,je,iz;
  double z,wj,az,swj,swjy,bii,thetai;
  /*double diff;*/

  for (i=1;i<=n;i++)
    {
      thetai=theta[i-1];
      ja=max(1,i-ih);
      je=min(n,i+ih);
      swj=0;
      swjy=0;


      for (j=ja;j<=je;j++)
	{
	  bii=bi[i-1]+bi[j-1];
	  z=thetai-theta[j-1];
	  z=z*z*bii*inv_lamakt;

	  if (z<100)
	    {
	      iz=z;
	      az=z-iz;
	      wj=(kerns[iz+1] - kerns[iz])*az + kerns[iz];
	      /*diff=i-j;*/
	      z=(double)(i-j)*inv_hakt;
	      z=z*z;

	      if (z<100)
		{	     
		  iz=z;
		  az=z-iz;
		  wj=wj*((kernl[iz+1] - kernl[iz])*az + kernl[iz]);
		  swj=swj+wj;
		  swjy = swjy+wj*y[j-1];
		}
	    }
	}
      ai[i-1]=swjy;
      bi[i-1]=swj;
    }
      
}



void gawsuni(double *y,
	     int *n,
	     double *hinit,
	     double *hincr,
	     double *hmax,
	     double *lamakt,
	     double *eta,
	     double *theta,
	     double *bi,
	     double *ai,
	     double *kernl,
	     double *kerns,
	     double *biold)
{
  int i;
  int ih;
  double z;
  double hakt;
  double inv_hakt;
  const int n_aux=*n;
  const double lamakt_aux=(2* *lamakt)/100;
  const double inv_lamakt=100/(2* *lamakt);
  const double hmax_aux=*hmax;
  const double eta_aux=*eta;
  const double hincr_aux=*hincr;
  const double hinit_aux=*hinit;

  hakt=hinit_aux*hincr_aux;

  inv_hakt=10/hakt;
  ih=hakt;
  lawsuni(y,n_aux,inv_hakt,ih,inv_lamakt,theta,bi,ai,kernl,kerns);
  i=n_aux;
  /*  for (i=0;i<n_aux;i++)*/
  while (i--)
    {
      z=eta_aux*(biold[i]*theta[i] - ai[i]) + ai[i];
      bi[i]=eta_aux * (biold[i] - bi[i]) + bi[i];
      theta[i]=z/bi[i];
    }
  memcpy(biold,bi,(n_aux*sizeof(bi[0])));
  hakt=hakt * hincr_aux;


  while (hakt<=hmax_aux)
    {
      inv_hakt=10/hakt;
      ih=hakt;
      lawsuni(y,n_aux,inv_hakt,ih,inv_lamakt,theta,bi,ai,kernl,kerns);
      i=n_aux;
      /*      for (i=0;i<n_aux;i++)*/
      while (i--)
	{
	  z=eta_aux*(biold[i]*theta[i] - ai[i]) + ai[i];
	  bi[i]=eta_aux * (biold[i] - bi[i]) + bi[i];
	  theta[i]=z/bi[i];
	}
      memcpy(biold,bi,(n_aux*sizeof(bi[0])));
      hakt=hakt * hincr_aux;
    }
}

 
