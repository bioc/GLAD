void iawsuni (double *y,
              int *n,
	      double *hinit,
	      double *bi,
	      double *ai,
	      double *kern,
	      double *theta);

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

