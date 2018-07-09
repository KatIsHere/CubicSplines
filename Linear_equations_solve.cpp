#include <cmath>

/* INPUT: A - array of pointers to rows of a square matrix having dimension N
*        Tol - small tolerance number to detect failure when the matrix is near degenerate
*  OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
*        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
*        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
*        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
*/
int LUPDecompose(double **A, const int& N, double Tol, int *P) {

	int i, j, k, imax;
	double maxA, *ptr, absA;

	for (i = 0; i <= N; i++)
		P[i] = i; //Unit permutation matrix, P[N] initialized with N

	for (i = 0; i < N; i++) {
		maxA = 0.0;
		imax = i;

		for (k = i; k < N; k++)
			if ((absA = fabs(A[k][i])) > maxA) {
				maxA = absA;
				imax = k;
			}

		if (maxA < Tol) throw "Matrix is degenerate, try another method"; //failure, matrix is degenerate

		if (imax != i) {
			//pivoting P
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			//pivoting rows of A
			ptr = A[i];
			A[i] = A[imax];
			A[imax] = ptr;

			//counting pivots starting from N (for determinant)
			P[N]++;
		}

		for (j = i + 1; j < N; j++) {
			A[j][i] /= A[i][i];

			for (k = i + 1; k < N; k++)
				A[j][k] -= A[j][i] * A[i][k];
		}
	}

	return 1;  //decomposition done 
}


/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
*  OUTPUT: x - solution vector of A*x=b
*/
void LUPSolve(double **A, int *P, double *b, const int& N, double *x) {

	for (int i = 0; i < N; i++) {
		x[i] = b[P[i]];

		for (int k = 0; k < i; k++)
			x[i] -= A[i][k] * x[k];
	}

	for (int i = N - 1; i >= 0; i--) {
		for (int k = i + 1; k < N; k++)
			x[i] -= A[i][k] * x[k];

		x[i] = x[i] / A[i][i];
	}
}

/* INPUT: A,P filled in LUPDecompose; N - dimension
*  OUTPUT: IA is the inverse of the initial matrix
*/
void LUPInvert(double **A, int *P, const int& N, double **IA) {

	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++) {
			if (P[i] == j)
				IA[i][j] = 1.0;
			else
				IA[i][j] = 0.0;

			for (int k = 0; k < i; k++)
				IA[i][j] -= A[i][k] * IA[k][j];
		}

		for (int i = N - 1; i >= 0; i--) {
			for (int k = i + 1; k < N; k++)
				IA[i][j] -= A[i][k] * IA[k][j];

			IA[i][j] = IA[i][j] / A[i][i];
		}
	}
}

/* INPUT: A,P filled in LUPDecompose; N - dimension.
*  OUTPUT: Function returns the determinant of the initial matrix
*/
double LUPDeterminant(double **A, int *P, const int& N) {

	double det = A[0][0];

	for (int i = 1; i < N; i++)
		det *= A[i][i];

	if ((P[N] - N) % 2 == 0)
		return det;
	else
		return -det;
}

double* BaireissSolve(double** A, double* b, const int & N) {
	double* x = new double[N];
	double** A_copy = new double*[N];
	int i, j, k, t;
	for (i = 0; i < N; ++i) {
		A_copy[i] = new double[N+1];
		for (j = 0; j < N; ++j) {
			A_copy[i][j] = A[i][j];
		}
		A_copy[i][N] = b[i];
	}
	
	double p0 = 1, p1 = A_copy[0][0];
	
	for (j = 0; j < N + 1; ++j) {
		for (i = 0; i < N; ++i) {
			if (i != j) {
				A_copy[i][j] = (A[i][j] * p1) / p0;
			}
		}
	}



	for (i = 0; i < N; ++i) {
		delete[]A_copy[i];
	}
	delete[]A_copy;

	return x;
}

double* GaussMethod(double** A, double*b, const int& n) {
	/* Finds a solution to the system of linear equations 
										using Gauss method*/
	double *x, max;
	int k, index;
	const double eps = 0.001;  // precision
	x = new double[n];
	k = 0;
	while (k < n)
	{
		// Find a row with max a[i][k]
		max = abs(A[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(A[i][k]) > max)
			{
				max = abs(A[i][k]);
				index = i;
			}
		}
		// Rearrange lines
		if (max < eps)
		{
			// no non-zero elements
			throw "Решение получить невозможно из-за нулевого столбца матрицы A";
		}
		for (int j = 0; j < n; j++)
		{
			double temp = A[k][j];
			A[k][j] = A[index][j];
			A[index][j] = temp;
		}
		double temp = b[k];
		b[k] = b[index];
		b[index] = temp;
		// Normalization of equations
		for (int i = k; i < n; i++)
		{
			double temp = A[i][k];
			if (abs(temp) < eps) continue; // skip for a zero element
			for (int j = 0; j < n; j++)
				A[i][j] = A[i][j] / temp;
			b[i] = b[i] / temp;
			if (i == k)  continue; // the equation does not deduct itself from itself
			for (int j = 0; j < n; j++)
				A[i][j] = A[i][j] - A[k][j];
			b[i] = b[i] - b[k];
		}
		k++;
	}
	// inverse substitution
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = b[k];
		for (int i = 0; i < k; i++)
			b[i] = b[i] - A[i][k] * x[k];
	}
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
