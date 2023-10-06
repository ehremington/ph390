#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*  Example  2.2                                                          *
%*  filename: ch02pr02.py                                                 *
%*  program listing number: 2.2                                           *
%*                                                                        *
%*  This program evaluates the derivative of a given function func(x)     *
%*  at x=1 using the mean finite difference method with the accuracy      *
%*  specified by tolerance.                                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/01/2017.                                    *
%**************************************************************************
"""

import numpy as np
import matplotlib.pyplot as plt
       
def func(x):  # define a function
    return x**3/3

if __name__ == "__main__":
    tol=input("Enter tolerance =")  # Read a tolerabce from the console
    tol=np.float(tol)
    x=1.0
    h=1.0   # initial interval
    diff_old=(func(x+h)-func(x-h))/(2.0*h)  #  derivative first try
    delta=np.finfo(float).max  # any value bigger than tol is OK.

    print("{0:^10} {1:^16} {2:^16}"
          .format('h','derivative','error'))
    while (delta>tol):
        h=h/2.0
        diff_new=(func(x+h)-func(x-h))/(2.0*h)  # improved derivative
        delta=np.abs(diff_new-diff_old)
        print("{0:10.3e} {1:16.10e} {2:16.10e}"
              .format(h,diff_new,delta))
        diff_old=diff_new
        
    print("Tolerance is OK.")

