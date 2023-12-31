#include <vector>
#include <map>

using namespace std;


extern "C"
{

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

  double IQRdiff(vector<double> vec);

  double IQR_vector_double(vector<double> vec);

  double quantile_vector_double(vector<double> vec, const double quantile);

  double median_fabs_double(const double *value, const int l);

  double mean_vector_double(vector<double> vec);

  double mad_vector_double(vector<double> vec);

  double var_vector_double(vector<double> vec, int unbiased=0);

  void printagg(vector<struct agg> agg_region);
  
  double kernelpen(double x, const double d);

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
			    const double *deletionValue,
			    const double *deltaNValue);


  void my_merge_int(const int *index_dest,
		    int *value_dest,
		    const int *index_src,
		    const int *value_src,
		    const int *length_dest,
		    const int *length_src);

  void detectOutliers(const double *LogRatio,
		      const int *Region,
		      int *OutliersAws,
		      int *OutliersMad,
		      int *OutliersTot,
		      const int *msize,
		      const double *alpha,
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
				      const double *deltaNValue,
				      //variables pour calcul la médiane par cluster
				      const double *LogRatio,
				      const int *NormalRange);


}
