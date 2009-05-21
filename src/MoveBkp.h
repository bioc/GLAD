extern "C"
{
  void MoveBkp_Delete_Bkp(const int *subBkpInfo_MoveBkp,
			  const int *subBkpInfo_PosOrder,
			  int *Breakpoints,
			  int *OutliersTot,
			  int *OutliersAws,
			  int *OutliersMad,
			  int *Level,
			  int *Region,
			  double *Smoothing,
			  int *GNL,
			  const int *l);

  void MoveBkp_updateOutliers (int *OutliersAws,
			       int *OutliersTot,
			       int *Level,
			       int *Region,
			       int *Breakpoints,
			       double *Smoothing,
			       int *ZoneGNL,
			       const int *l);




  void MoveBkp_Step1(const int *subBkpInfo_MoveBkp,
		     const int *subBkpInfo_PosOrder,
		     const double *LogRatio,
		     double *NextLogRatio,
		     const int *Chromosome,
		     const int *PosOrder,
		     int *Breakpoints,
		     int *OutliersTot,
		     int *OutliersAws,
		     int *OutliersMad,
		     int *Level,
		     int *Region,
		     double *Smoothing,
		     int *GNL,
		     int *NormalRange,
		     const double *NormalRef,
		     const double *deltaN,
		     const int *lensubBkp,
		     const int *l);

  void MoveBkp_Step2(int *OutliersAws,
		     int *OutliersTot,
		     int *Level,
		     int *Region,
		     int *Breakpoints,
		     // variables pour faire la jointure
		     int *ZoneGNL,
		     int *value_dest,
		     const int *length_dest,
		     double *Smoothing,
		     const double *forceGL1Value,
		     const double *forceGL2Value,
		     const double *NormalRefValue,
		     const double *ampliconValue,
		     const double *deletionValue,
		     const double *deltaN,
		     //variables pour calcul la médiane par cluster
		     const double *LogRatio,
		     int *NormalRange);

  void rangeGainLoss(const double *Smoothing,
		     const int *ZoneGNL,
		     const int *OutliersTot,
		     double *minG,
		     double *maxL,
		     double *minAmp,
		     double *maxDel,
		     const int *l);

  int  testSingle(const double LogRatio,
		  const double NextLogRatio,
		  const double Smoothing,
		  const double SmoothingNext);


}
