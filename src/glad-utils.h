#include <vector>
#include <map>

using namespace std;


extern "C"
{

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
			    const double *deletionValue);


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



}
