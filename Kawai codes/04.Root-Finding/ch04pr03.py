#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  4.3                                                       *
%*     filename: ch04pr03.py                                               *
%*     program listing number: 4.3                                        *
%*                                                                        *
%*     This program finds roots of a cubic equation                       *
%*               a*x^3 + b*x^2 + c*x + d = 0                              *
%*     using the bisection method.                                        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/24/2018.                                    *
%**************************************************************************
"""
import numpy as np

# define a cubic equation
def f(x):
    return x**3-9.0*x**2+23.0*x-15

if __name__ == "__main__": 
    
    tol=input("Enter tolerance =")  # Read a tolerabce from the console
    epsilon=np.float(tol)
    
    # initial bracket
    x1=1.5
    x2=4.0
    f1 = f(x1)
    f2 = f(x2)

    # iteration counter
    n=0

    # mid point
    xm = (x1+x2)/2.0
    fm = f(xm)
    
    while x2-x1 > epsilon:
        if f1*fm < 0:  # root in the lower half
            x2=xm
            f2=fm
        else:          # root in the upper half
            x1=xm
            f1=fm

        xm = (x1+x2)/2.0      # new mid point
        fm = f(xm)
        n+=1

    print("Answer = {0:7.5f}, (iteration = {1:3d})"
          .format(xm,n))


