# CubicSplines
A program for calculating natural smoothing cubic splines - an effective tool in function interpolation. Written in C++

----------------------------------------------------------------------
Materials about splines can be found on wikipedia: https://en.wikipedia.org/wiki/Spline_(mathematics)

Basically, it interpolates function f(x) not on the whole [a, b] section, but on subsection [X[i], X[i+1]], 
for a = X[0] <= X[1] <= .... <= X[n] = b

This version can also create smoothing cubic splines, which help to interpolate functions with a lot of noises.
