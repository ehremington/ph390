#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  8.6                                                       *
%*     filename: ch08pr06.py                                              *
%*     program listing number: 8.6                                        *
%*                                                                        *
%*     This program solves a simple linear equation with LU decomposition.*
%*     MATLAB function lu() is used.                                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2017.                                    *
%**************************************************************************
"""

import numpy as np
import scipy.linalg as la

# Define a matrix
A=np.matrix([[3., -1., 4.],[2., 0., -1.],[0., 3., 2.]])
b=np.matrix([2.,-1.,3.]).transpose()
x=np.matrix(np.zeros(3)).transpose()
y=np.matrix(np.zeros(3)).transpose()

# LU dcomposition
P, L, U = la.lu(A)
P=np.matrix(P)
U=np.matrix(U)
L=np.matrix(L)

# Rcover the original matrix
S=P*L*U
# Show the results
print('\nA (Original Matrix)')
print(A)
print('\nL (Lower Triangular Matrix)')
print(L)
print('\nU (Upper Triangular Matrix)')
print(U)
print('\nP (Permutation Matrix)')
print(P)
print('\nP*L*U')
print(S)

b = P*b
# forward substition
for i in range(0,3):
    Ly=0.0
    for j in range(0,i):
        Ly = Ly+L[i,j]*y[j]

    y[i] = (b[i]-Ly)/L[i,i]

# backsubstitution
for i in range(2,-1,-1):
    Ux=0.0
    for j in range(i+1,3):
        Ux = Ux+U[i,j]*x[j]

    x[i] = (y[i]-Ux)/U[i,i]

print('\nSolution x')
print(x)




