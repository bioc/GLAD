#include <vector>
#include <map>

using namespace std;

extern "C"
{

  void OptmisationBreakpointsStep(const int *Chromosome,
				  double *Smoothing,
				  int *NormalRange,
				  const double *NormalRef,
				  const double *deltaN,
				  // variable pour loop_chromosome_removeLevel
				  const double *LogRatio,
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
				  const int *BkpDetected,
				  const int *l); // nombre total de sondes

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

  void updateBkpRL (int *Region,
		    int *OutliersAws,
		    int *Breakpoints,
		    const int *PosOrder,
		    double *NextLogRatio,
		    const double *LogRatio,
		    const int *l);

  void makeRegionLevelID(const int *Chromosome,
			 int *Level,
			 const int n);

  double computeLike(vector<struct agg> agg_region, double lambda, double sumkernelpen);

  double computeSumKernelPen(vector<struct agg> agg_region, double sigma, double d);


}
