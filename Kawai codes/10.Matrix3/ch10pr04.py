#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 10.6                                                       *
%*     filename: ch10pr04.py                                              *
%*     program listing number: 10.4                                       *
%*                                                                        *
%*     This program finds eigenvalues and eigenvectos of a symmetric      *
%*     matrix using the QR algorithm method.                              *
%*                                                                        *
%*     Uses NUMPY function qr()                                           *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/18/2017.                                    *
%**************************************************************************
"""
import numpy as np

# Define the matrix
A=np.matrix([[1.,-4.,2.],[-4.,1.,-2.],[2.,-2.,-2.]])

# Tolerance
tol = 1e-8;

# Magnitude of total off-diagonal elemens
S=0.0
for i in range(0,3):
    for j in range(i+1,3):
        S=S+A[i,j]**2

# Initial transformation matrix
P=np.matrix(np.identity(3))

n=0
while S > tol:
    n+=1
    [Q, R] =np.linalg.qr(A)  # QR decomposition
    A=R*Q  # Orthogonal transformation
    P=P*Q  # Accumulating the transformation
    # Error evaluation
    S=0.0
    for i in range(0,3):
        for j in range(i+1,3):
            S=S+A[i,j]**2

print('# of iterations ={0:d}'.format(n))

print('\nTransformed Matrix\n')
print(A)

print('\nTransformation Matrix\n')
print(P)
