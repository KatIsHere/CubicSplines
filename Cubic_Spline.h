#pragma once
double* Smoothing_Spine(double* X, double* F_wrong, double* P, const int& N);
double Cubic_Spline(double x, double* X, double* F, double* M, const int& N);
double* Cubic_Spline_Vector(double** X, double** F, const int& N);