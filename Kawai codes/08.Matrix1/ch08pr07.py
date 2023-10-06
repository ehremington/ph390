#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  8.7                                                       *
%*     filename: ch08pr07.py                                              *
%*     program listing number: 8.7                                        *
%*                                                                        *
%*     This program solves a tridiagonal system with backward elimination.*
%*     Then, find solution by forward substitution.                       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2017.                                    *
%**************************************************************************
"""
import numpy as np

# Define matrices. No need to use the full matrix.
d=[2,3,4,3] # diagonal elements
u=[2,3,3,0] # above diagonal
l=[0,2,3,3] # below diagonal
b=[1,2,3,4] # right hand side
D = np.zeros(4)
Y = np.zeros(4)
Z = np.zeros(4)
s = np.zeros(4)
x = np.zeros(4)
# Calculation of determinant
D[0]=d[0]
D[1]=d[1]*D[0]-l[0]*u[0]
for i in range(2,4):
    D[i]=d[i]*D[i-1]-l[i-1]*u[i-1]*D[i-2]


print('Determinant {0:8.3f}='.format(D[3]))
if D[3] == 0:
    exit('Singular')

# Decomposition by backword elimination
Y[2]=-l[3]/d[3]
Z[2]= b[3]/d[3]
for i in range(2,0,-1):
    Y[i-1]=-l[i]/(d[i]+u[i]*Y[i])
    Z[i-1]=(b[i]-u[i]*Z[i])/(d[i]+u[i]*Y[i])


# Forward substitution
x[0]=(b[0]-u[0]*Z[0])/(d[0]+u[0]*Y[0])
for i in range(0,3):
    x[i+1]=Y[i]*x[i]+Z[i]

# Answer
print('Solution x')
print(x)

# Check the errors. 
s[0]=          d[0]*x[0]+u[0]*x[1]-b[0]
s[1]=l[1]*x[0]+d[1]*x[1]+u[1]*x[2]-b[1]
s[2]=l[2]*x[1]+d[2]*x[2]+u[2]*x[3]-b[2];
s[3]=l[3]*x[2]+d[3]*x[3]-b[3]
print('Error')
print(s)
