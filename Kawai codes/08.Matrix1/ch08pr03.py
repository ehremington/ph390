#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  8.3                                                       *
%*     filename: ch08pr03.py                                              *
%*     program listing number: 8.3                                        *
%*                                                                        *
%*     This program solves a simple linear equation with the Gaussian     *
%*     eliminationa and backsubstitution methods.                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/06/2017.                                    *
%**************************************************************************
"""
import numpy as np

# Set a linear equation
N=3;
A=np.matrix([[3.,-1.,4.],[2.,0.,-1.],[0.,3.,2.]])
b=np.matrix.transpose(np.matrix([2.,-1.,3.]))
x=np.matrix.transpose(np.matrix(np.zeros(N)))

#forward elimination
for n in range(0,N-1):
  
    for i in range(n+1,N):
        M=-A[i,n]/A[n,n]
        A[i,n+1:N]=M*A[n,n+1:N]+A[i,n+1:N]
        b[i]=M*b[n]+b[i]

    A[n+1,n]=0.0 


# backsubstitution
for i in range(2,-1,-1):
    Ax=0.0
    for j in range(i+1,3):
        Ax = Ax+A[i,j]*x[j]

    x[i] = (b[i]-Ax)/A[i,i]


# result
print('\nA=\n')
print(A)
print('\nb=\n')
print(b)
print('\nx=\n')
print(x)