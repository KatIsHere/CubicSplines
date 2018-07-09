#pragma once
double* GaussMethod(double** A,double* b, const int& m);
double* TridiagonalSolve(const double *a, const double *b, double *c, double *d, const int& n);
double* TridiagonalSolve(const int *a, const int *b, double *c, int *d, const int& n);
int LUPDecompose(double **A, const int& N, double Tol, int *P);
void LUPSolve(double **A, int *P, double *b, const int& N, double *x);
void LUPInvert(double **A, int *P, const int& N, double **IA);
double LUPDeterminant(double **A, int *P, const int& N);
