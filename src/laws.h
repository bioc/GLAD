/* int max0(int N1, int N2); */

/* int min0(int N1, int N2); */


/* void lawsuni(double *y, */
/* 	     const int n, */
/* 	     double inv_hakt, */
/* 	     int ih, */
/* 	     const double lamakt, */
/* 	     double *theta, */
/* 	     double *bi, */
/* 	     double *ai, */
/* 	     double *kernl, */
/* 	     double *kerns); */



/* struct param */
/* { */
/*   double *y; */
/*   int n; */
/*   double *kerns; */
/*   double *kernl; */
/*   double *ai; */
/*   double *bi; */
/*   double *theta; */
/* }; */


/* void lawsuni(struct param *par, */
/* 	     double inv_hakt, */
/* 	     int ih, */
/* 	     const double inv_lamakt); */

void lawsuni(double *y,
	     const int n,
	     double inv_hakt,
	     int ih,
	     const double inv_lamakt,
	     double *theta,
	     double *bi,
	     double *ai,
	     double *kernl,
	     double *kerns);



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
	     double *biold);
