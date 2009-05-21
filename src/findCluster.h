#include <vector>
#include <map>


using namespace std;

extern "C"
{
  int clusterglad(map<int, struct agg> map_clusterRegion,
		  int *treemerge,
		  const int min,
		  const int max,
		  const double sigma,
		  const double d,
		  const double lambda);


  void mergeLike(map<int, struct agg> map_clusterRegion,
		 double *logVar,
		 double *Mean,
		 const int *classe,
		 const int which_classe);

  void findCluster(const double *LogRatio,
		   const int *Region,
		   const int *OutliersTot,
		   int *zone,
		   int *method,
		   // paramètres pour clusterglad
		   const double *sigma,
		   const double *d,
		   const double *lambda,
		   const int *nmin,
		   const int *nmax,
		   int *nbclasses,
		   const int *l);
}
