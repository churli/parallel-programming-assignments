
void compute(unsigned long **a, unsigned long **b, unsigned long **c, unsigned long **d, int N, int num_threads) {

	// perform loop fusion to transform this loop and parallelize it with OpenMP
	#pragma omp parallel for num_threads(num_threads)
	for (int i = 1; i < N; i++) {
		a[i][1] = 2 * b[i][1];
		d[i][1] = a[i][1] * c[i][1];
		int j;
		for (j = 2; j < N; j++) {
			a[i][j] = 2 * b[i][j];
			d[i][j] = a[i][j] * c[i][j];
			c[i][j - 1 - 1] = a[i][j - 1 - 1] - a[i][j + 1 - 1];
		}
		c[i][j - 1 - 1] = a[i][j - 1 - 1] - a[i][j + 1 - 1];
	}
}
