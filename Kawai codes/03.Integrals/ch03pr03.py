#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*  Example  3.4                                                          *
%*  filename: ch03pr03.py                                                 *
%*  program listing number: 3.3                                           *
%*                                                                        *
%*  This program numerically integrates x^3*exp(x)/(exp(x)-1) from        *
%*  x=0 to infinity using 8-point Gaussian Laguerre Quadrature.           *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/12/2017.                                    *
%**************************************************************************
"""
import numpy as np

def f(x):
    return x**3*np.exp(x)/(np.exp(x)-1.)

if __name__ == "__main__":    
    N=8
# evaluation points and weights for 8 point Gaussian quadrature
    x=np.array([1.7027963230510100e-1, 9.0370177679937991e-1,
                2.2510866298661307,    4.2667001702876588,
                7.0459054023934657,    1.0758516010180995e+1,
                1.5740678641278005e+1, 2.2863131736889264e+1])

    w=np.array([3.6918858934163753e-1, 4.1878678081434296e-1,
                1.7579498663717181e-1, 3.3343492261215652e-2,
                2.7945362352256725e-3, 9.0765087733582131e-5,
                8.4857467162725315e-7, 1.0480011748715104e-9])

    gauss=(w*f(x)).sum()  #Gaussian quadrature
    exact=np.pi**4/15.
    print("{0:3d} point Gaussian Laguerre Quadrature".format(N))
    print("  Exact={0:18.12e}\n  Gauss={1:18.12e}\n  Error={2:18.12e}"
          .format(exact,  gauss, abs(exact-gauss)))
