double* GaussMethod(double** A, double* b, const int& m) {
	/* Finds a solution to the system of linear equations 
										using Gauss method*/
	// Ax = b, 
	// A: m x m, b: m

	double* x = new double[m];
	double** A_copy = new double*[m];
	double* b_copy = new double[m];
	int i, j, k, t;
	for (i = 0; i < m; ++i) {
		A_copy[i] = new double[m];
		b_copy[i] = b[i];
		for (j = 0; j < m; ++j) {
			A_copy[i][j] = A[i][j];
		}
	}
	double ai;
	for (i = 0; i < m; ++i) {
		ai = A_copy[i][i];
		if (ai == 0) continue;
		for (j = 0; j < m; ++j) {
			A_copy[i][j] /= ai;
		}
		b_copy[i] /= ai;
		for (k = m - 1; k > i; --k) {
			for (j = 0; j < m; ++j) {
				A_copy[k][j] -= (A_copy[i][j] * A_copy[k][i]);
			}
			b_copy[k] -= (b_copy[i] * A_copy[k][i]);
		}
	}

	x[m - 1] = b_copy[m - 1];
	double suma;
	for (i = m - 2; i >= 0; --i) {
		suma = 0.0;
		for (j = m - 1; j > i; --j) {
			suma += A_copy[i][j] * x[j];
		}
		x[i] = b_copy[i] - suma;
	}

	for (i = 0; i < m; ++i) {
		delete[]A[i];
	}
	delete[]A_copy;
	delete[]b_copy;

	return x;
}

double* GaussMethod(int** A, int* b, const int& m) {
	/* Finds a solution to the system of linear equations
	using Gauss method*/
	// Ax = b, 
	// A: m x m, b: m

	double* x = new double[m];
	double** A_copy = new double*[m];
	double* b_copy = new double[m];
	int i, j, k, t;
	for (i = 0; i < m; ++i) {
		A_copy[i] = new double[m];
		b_copy[i] = b[i];
		for (j = 0; j < m; ++j) {
			A_copy[i][j] = (double)A[i][j];
		}
	}
	double ai;
	for (i = 0; i < m; ++i) {
		ai = A_copy[i][i];
		if (ai == 0) continue;
		for (j = 0; j < m; ++j) {
			A_copy[i][j] /= ai;
		}
		b_copy[i] /= ai;
		for (k = m - 1; k > i; --k) {
			for (j = 0; j < m; ++j) {
				A_copy[k][j] -= (A_copy[i][j] * A_copy[k][i]);
			}
			b_copy[k] -= (b_copy[i] * A_copy[k][i]);
		}
	}

	x[m - 1] = b_copy[m - 1];
	double suma;
	for (i = m - 2; i >= 0; --i) {
		suma = 0.0;
		for (j = m - 1; j > i; --j) {
			suma += A_copy[i][j] * x[j];
		}
		x[i] = b_copy[i] - suma;
	}

	for (i = 0; i < m; ++i) {
		delete[]A[i];
	}
	delete[]A_copy;
	delete[]b_copy;

	return x;
}

double* TridiagonalSolve(const double *a, const double *b, double *c, double *d, const int& n) {
	// a[0] == 0
	// c[n-1] == 0
	double* x = new double[n];
	/* Modify the coefficients. */
	c[0] /= b[0];	
	d[0] /= b[0];	
	double id;	int i;
	for (i = 1; i < n; i++) {
		id = 1 / (b[i] - c[i - 1] * a[i]);  /* Division by zero risk. */
		c[i] *= id;	                         /* Last value calculated is redundant. */
		d[i] = (d[i] - d[i - 1] * a[i]) * id;
	}

	x[n - 1] = d[n - 1];
	for (i = n - 2; i >= 0; i--)
		x[i] = d[i] - c[i] * x[i + 1];

	return x;
}

double* TridiagonalSolve(const int *a, const int *b, double *c, int *d, const int& n) {
	// a[0] == 0
	// c[n-1] == 0
	double* x = new double[n];
	/* Modify the coefficients. */
	c[0] /= b[0];
	d[0] /= b[0];
	double id;	int i;
	for (i = 1; i < n; i++) {
		id = 1 / (b[i] - c[i - 1] * a[i]);  /* Division by zero risk. */
		c[i] *= id;	                         /* Last value calculated is redundant. */
		d[i] = (d[i] - d[i - 1] * a[i]) * id;
	}

	x[n - 1] = d[n - 1];
	for (i = n - 2; i >= 0; i--)
		x[i] = d[i] - c[i] * x[i + 1];

	return x;
}