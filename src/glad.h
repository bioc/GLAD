#include <vector>
#include <map>

extern "C"
{
  using namespace std;

  double IQRdiff(vector<double> vec);

  double IQR_vector_double(vector<double> vec);

  double quantile_vector_double(vector<double> vec, const double quantile);

  double median_fabs_double(const double *value, const int l);

  double mean_vector_double(vector<double> vec);

  double mad_vector_double(vector<double> vec);

  double var_vector_double(vector<double> vec, int unbiased=0);

  double computeSumKernelPen(vector<struct agg> agg_region, double sigma, double d);

  double computeLike(vector<struct agg> agg_region, double lambda, double sumkernelpen);

  void printagg(vector<struct agg> agg_region);
  
  double kernelpen(double x, const double d);

  void detectOutliers(const double *LogRatio,
		      const int *Region,
		      int *OutliersAws,
		      int *OutliersMad,
		      int *OutliersTot,
		      const int *msize,
		      const double *alpha,
		      const int *l);

  void my_merge(const int *index_dest,
		double *value_dest,
		const int *index_src,
		const double *value_src,
		const int *length_dest,
		const int *length_src);

  void my_merge_int_forceGL(const int *index_dest,
			    int *value_dest,
			    const int *index_src,
			    const int *value_src,
			    const int *length_dest,
			    const int *length_src,
			    const double *Smoothing,
			    const double *forceGL1Value,
			    const double *forceGL2Value,
			    const double *NormalRefValue,
			    const double *ampliconValue,
			    const double *deletionValue);


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

  void putLevel(double *Smoothing,
		const double *LogRatio,
		int *Level,
		int *nblevel,
		const int *l);


}
