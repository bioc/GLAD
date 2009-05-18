extern "C"
{

  void awsBkp(const double *Smoothing,
	      int *OutliersAws,
	      int *Level,
	      int *nbregion,
	      int *regionChr,
	      int *Breakpoints,
	      int *bkp_detected,
	      const int *l);

  void putLevel_awsBkp(//variables pour putLevel
		       double *Smoothing,
		       const double *LogRatio,
		       int *Level,
		       int *nblevel,
		       const int *l,
		       // variables pour awsBkp
		       int *OutliersAws,
		       int *nbregion,
		       int *regionChr,
		       int *Breakpoints,
		       int *bkp_detected);

  void putLevel(double *Smoothing,
		const double *LogRatio,
		int *Level,
		int *nblevel,
		const int *l);

}
