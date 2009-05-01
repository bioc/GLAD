#include <vector>
#include <map>

extern "C"
{
  using namespace std;

  double median_vector_double(vector<double> vec);

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




}
