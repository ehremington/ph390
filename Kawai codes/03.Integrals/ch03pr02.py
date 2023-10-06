#!/usr/bin/env python3
"""
%**************************************************************************
%*  Example  3.3                                                          *
%*  filename: ch03pr02.py                                                 *
%*  program listing number: 3.2                                           *
%*                                                                        *
%*  This program integrates 1/(sqrt(x)*(1+x)) from x=0 to x=1             *
%*  by removing singularity at x=0. Trapezoidal rule is used              *
%*  for the proper part of integral.                                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/11/2017.                                    *
%**************************************************************************
"""
import numpy as np

a = 1.0 # upper bound
N = 100 # number of segments
h = a/N # width of segments
# integration of sqrt(x)/(1+x) with trapezoidal rule
S = np.sqrt(a)/(1.0+a)/2.0; # bundary value devided by 2
for i in range(1,N):
    x = i*h
    f = np.sqrt(x)/(1+x)
    S = S +f

proper = S*h # integral of proper part
singular = 2*np.sqrt(a)
# singular part
total = singular - proper
exact = np.pi - 2*np.arctan(1/np.sqrt(a))
print("Numerical = {0:18.12e}".format(total))
print("    Exact = {0:18.12e}".format(exact))

