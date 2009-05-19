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
		  const int *l);

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
				   const int *BkpDetected);

  double computeLike(vector<struct agg> agg_region, double lambda, double sumkernelpen);

  double computeSumKernelPen(vector<struct agg> agg_region, double sigma, double d);


}
