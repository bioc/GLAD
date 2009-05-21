extern "C"
{
  void updateFilterBkp(const int *Chromosome,
		       int *Breakpoints,
		       int *Level,
		       const int *PosOrder,
		       double *NextLogRatio,
		       const double *LogRatio,
		       const int *maxLevel,
		       // ajout de variables pour updateOutliers
		       int *OutliersAws,
		       double *Smoothing,
		       // ajout de variables pour detectOutliers
		       int *OutliersMad,
		       int *OutliersTot,
		       const int *msize,
		       const double *alpha,
		       const int *l,
		       const double *NormalRef,
		       const double *deltaN,
		       int *NormalRange);    
}
