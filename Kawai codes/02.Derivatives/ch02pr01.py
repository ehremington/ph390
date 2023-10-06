#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*  Example  2.1                                                          *
%*  filename: ch02pr01.py                                                 *
%*  program listing number: 2.1                                           *
%*                                                                        *
%*  This program evaluates the derivative of a given function func(x)     *
%*  at x=1 using the three finite difference methods.                     *
%*  Errors in forward, backward and mean value methods are plotted.       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  12/31/2016.                                    *
%**************************************************************************
"""

import numpy as np
import matplotlib.pyplot as plt
       
def func(x):
    return x**3/3

if __name__ == "__main__":
    x=1.0
    imax=20
    h=np.zeros(imax)
    d_f=np.zeros(imax)
    d_b=np.zeros(imax)
    d_m=np.zeros(imax)
    err_f=np.zeros(imax)
    err_b=np.zeros(imax)
    err_m=np.zeros(imax)
    i=0
    print("{0:^62}".format('Absolute Errors'))
    print("{0:^6} {1:^18} {2:^18} {3:^20}"
          .format('h','forward','backward','mean value'))
    while(i<imax):
    # Small displacement
        h[i]=10**(-i)
    
    # Evaluation of numerical derivative
        d_f[i]=(func(x+h[i])-func(x))/h[i] # Forward diffrence
        d_b[i]=(func(x)-func(x-h[i]))/h[i] # Backward diffrence
        d_m[i]=(func(x+h[i])-func(x-h[i]))/(2*h[i]) # Mean value

    # Errors
        err_f[i]=abs(1.-d_f[i])
        err_b[i]=abs(1.-d_b[i])
        err_m[i]=abs(1.-d_m[i])
        print("{0:6.1e} {1:18.10e} {2:18.10e} {3:18.10e}"
              .format(h[i],err_f[i],err_b[i],err_m[i]))
        i=i+1

    # Plot data
    plt.ioff()
    plt.figure(figsize=(12,5))
    plt.subplot(1,2,1)
    plt.semilogx(h,d_f, '--ob', label='forward')
    plt.semilogx(h,d_b, '--dg', label='backword')
    plt.semilogx(h,d_m, '--sr', label='mean')
    plt.legend(loc=2)
    plt.xlabel('h')
    plt.ylabel('Derivative')
    
    plt.subplot(1,2,2)
    plt.loglog(h,err_f, '--ob', label='forward')
    plt.loglog(h,err_b, '--dg', label='backword')
    plt.loglog(h,err_m, '--sr', label='mean')
    plt.legend(loc=3)
    plt.xlabel('h')
    plt.ylabel('Absolute error')
    plt.show()