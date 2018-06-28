#pragma once
double* GaussMethod(double** A, double* b, const int& m);
double* GaussMethod(int** A, int* b, const int& m);
double* TridiagonalSolve(const double *a, const double *b, double *c, double *d, const int& n);
double* TridiagonalSolve(const int *a, const int *b, double *c, int *d, const int& n);