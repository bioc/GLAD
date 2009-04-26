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

}
