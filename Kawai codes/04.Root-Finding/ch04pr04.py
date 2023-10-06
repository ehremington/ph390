#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 4.4                                                        *
%*     filename: ch04pr04.m                                               *
%*     program listing number: 4.4                                        *
%*                                                                        *
%*     This program seeks the root of                                     *
%*             cos(3*x)*sin(x)=0                                          *
%*     between x=0.2 and 0.8 using bisection, Newton-Raphson, and         *
%*     secant methods.                                                    *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/14/2017.                                             *
%**************************************************************************
"""
import sys
import numpy as np

# define function and its derivative
def f(x):
    return np.cos(3.0*x)*np.sin(x)

def df(x):
    return -3.0*np.sin(3*x)*np.sin(x) + np.cos(3.0*x)*np.cos(x);

if __name__ == "__main__": 
    # initial bracket
    x1, x2 = input("Initial blacket (a pair of numbers) = ").split()
    x1, x2 = [np.float(x1), np.float(x2)]
    
    f1=f(x1)
    f2=f(x2)
    if f1*f2 > 0:
        sys.exit('Bracket is incorrect')
        
    # tolerance
    epsilon = input("Tolerance = ")
    epsilon = np.float(epsilon)

    # bisection 10 iterations
    n=0 # iteration counter

    xm = (x1+x2)/2.0
    fm = f(xm)
    
    while n<10:
        if f1*fm < 0:  # root in the lower half
            x2=xm
            f2=fm
        else:          # root in the upper half
            x1=xm
            f1=fm

        xm = (x1+x2)/2.0  # new mid point
        fm = f(xm)
        n+=1

    print("Bisection      = {0:12.8f}  (iteration= {1:4d})"
          .format(xm,n))

    # Newton-Raphson method
    x = xm
    fx = fm
    n = 0
    while abs(fx)> epsilon:
        dfx = df(x)
        x = x - fx/dfx
        fx = f(x)
        n+=1

    print("Newton-Raphson = {0:12.8f}  (iteration=  {1:4d})"
          .format(x,n))

    # Secant method
    dx = (x2-x1)/10.
    x1 = xm
    f1 = fm
    x2 = x1 + dx
    f2 = f(x2)
    n = 0
    while abs(f2)> epsilon:
        x = x2 - (x2-x1)/(f2-f1)*f2
        x1 = x2
        f1 = f2
        x2 = x
        f2 = f(x)
        n+=1

    print("Secant         = {0:12.8f}  (iteration=  {1:4d})"
          .format(x,n))
