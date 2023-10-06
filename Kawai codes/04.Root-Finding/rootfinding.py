#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*   Secion 4.6.                                                          *
%*   filename: rootfinding.py                                             *
%*   program listing number: 4.6-2                                        *
%*                                                                        *
%*     Inputs                                                             *
%*         f = function name                                              *
%*         x1, x2 = bracket                                               *
%*         N = interations of bisection                                   *
%*         tol = tolerance for scant method                               *
%*     Output:                                                            *
%*         x = root                                                       *
%*                                                                        *
%*     This program find a root of a given funcion f(x)=0                 *
%*     using the bisection and secant root finding method                 *
%*     finding methods.                                                   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
"""
def findroot(f,x1,x2,N,tol):
    
    f1=f(x1)
    f2=f(x2)
    if f1*f2 > 0:
        exit('Bracket is incorrect')

    # Bisection method (N iteration)
    n=0  # iteration counter

    # mid point
    xm = (x1+x2)/2.0
    fm = f(xm)

    while n<N+1:
        if f1*fm < 0.0:  # root in the lower half
            x2=xm
            f2=fm
        else:            # root in the upper half
            x1=xm
            f1=fm

        xm = (x1+x2)/2.0  # new mid point
        fm = f(xm)
        n+=1

    # Secant method
    dx = (x2-x1)/10.0
    x1 = xm
    f1 = fm
    x2 = x1 + dx
    f2 = f(x2)
    n = 0
    while abs(f2)> tol:
        x = x2 - (x2-x1)/(f2-f1)*f2
        x1 = x2
        f1 = f2
        x2 = x
        f2 = f(x)
        n+=1

    return x