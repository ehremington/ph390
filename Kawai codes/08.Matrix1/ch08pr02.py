#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  8.2                                                       *
%*     filename: ch08pr02.py                                              *
%*     program listing number: 8.2                                        *
%*                                                                        *
%*     This program solves a upper-triangular linear equation with        *
%*     the backsubstitution method.                                       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/06/2017.                                    *
%**************************************************************************
"""
import numpy as np

# define matrix A and vector b
A=np.matrix([[3, -1, 4],[0,2,-1],[0,0,2]])
b=np.matrix([-1,-2,4]).transpose()
x=np.matrix([-1,-2,4]).transpose()

# backsubstitution
for i in range(2,-1,-1):
    Ax=0.0
    for j in range(i+1,3):
        Ax = Ax+A[i,j]*x[j]

    x[i] = (b[i]-Ax)/A[i,i]

print(x)

