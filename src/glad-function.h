extern "C"
{

  int  testSingle(const double LogRatio,
		  const double NextLogRatio,
		  const double Smoothing,
		  const double SmoothingNext);

  void quicksort_int(int* array, int *order, int left, int right);

  int partition(int* array, int *order, int left, int right);

  int findMedianOfMedians(int* array, int *order, int left, int right);

  int findMedianIndex(int* array, int *order, int left, int right, int shift);

  void swap(int * array, int *order, int a, int b);


  void updateLevel (const int *Chromosome,
		    const int *Breakpoints,
		    int *Level,
		    const int *PosOrder,
		    double *NextLogRatio,
		    const double *LogRatio,
		    const int *maxLevel,
		    const int *l);


  void updateOutliers (int *OutliersAws,
                       int *Level,
                       int *Breakpoints,
		       double *Smoothing,
                       const int *l);

  void updateBkpRL (int *Region,
		    int *OutliersAws,
		    int *Breakpoints,
		    const int *PosOrder,
		    double *NextLogRatio,
		    const double *LogRatio,
		    const int *l);

  void compute_median_smoothing(const double *LogRatio,
				const int *ByValue,
				double *Smoothing,
				const int *l);

  void compute_NormalRange(const double *Smoothing,
			   const double *NormalRef,
			   const int *Level,
			   int *NormalRange,
			   const double *deltaN,
			   const int *l);



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

  void compute_cluster_LossNormalGain(// variables pour faire la jointure
				      const int *ZoneGen,
				      int *value_dest,
				      const int *length_dest,
				      const double *Smoothing,
				      const double *forceGL1Value,
				      const double *forceGL2Value,
				      const double *NormalRefValue,
				      const double *ampliconValue,
				      const double *deletionValue,
				      //variables pour calcul la médiane par cluster
				      const double *LogRatio,
				      const int *NormalRange);



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


}
