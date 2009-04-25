double kernelpen(double value, const double d);

int  testSingle(const double LogRatio,
		const double NextLogRatio,
		const double Smoothing,
		const double SmoothingNext);

void quicksort_int(int* array, int *order, int left, int right);

int partition(int* array, int *order, int left, int right);

int findMedianOfMedians(int* array, int *order, int left, int right);

int findMedianIndex(int* array, int *order, int left, int right, int shift);

void swap(int * array, int *order, int a, int b);
