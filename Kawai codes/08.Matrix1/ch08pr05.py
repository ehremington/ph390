#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  8.5                                                       *
%*     filename: ch08pr05.py                                              *
%*     program listing number: 8.5                                        *
%*                                                                        *
%*     This program calculates the inverse of a given matrix using        *
%*     Gaussi-Jordin methods.                                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2017.                                    *
%**************************************************************************
"""
import numpy as np

# Set a linear equation
N=3
A0=np.matrix([[3.,-1.,4.],[2.,0.,-1.],[0.,3.,2.]])
A=np.identity(N)
A[:,:]=A0[:,:]
b=np.matrix(np.identity(N))  # permutation matrix is initially an identity matrix
x=np.matrix(np.identity(N))
S=np.zeros(N)
TMP2=np.zeros(N)

# scale factors
for i in range(0,N): 
    S[i]=A[i,:].max()


for n in range(0,N-1):
    # Look for the pivot row
    j=n
    Amax=abs(A[n,n]/S[n])
    for i in range(n,N):
        AS=abs(A[i,n]/S[i])
        if AS > Amax:
            j=i
            Amax = AS

    # Carry out pivoting
    if j != n:
        for i in range(n,N): 
            TMP=A[n,i]
            A[n,i]=A[j,i]
            A[j,i]=TMP

        TMP2[:]=b[n,:]
        b[n,:]=b[j,:]
        b[j,:]=TMP2[:]

    # Gaussian elimination
    for i in range(n+1,N): 
        M=-A[i,n]/A[n,n]
        A[i,n+1:N]=M*A[n,n+1:N]+A[i,n+1:N]
        b[i,:]=M*b[n,:]+b[i,:]

    A[n+1,n]=0.0  

# backsubstitution
Ax=np.zeros(N)
for i in range(N-1,-1,-1):
    Ax=np.zeros(N)
    for j in range(i+1,N):
        for k in range(0,N):
            Ax[k] = Ax[k]+A[i,j]*x[j,k]

    x[i,:] = (b[i,:]-Ax[:])/A[i,i]

# result
print('\nInvers of A=')
print(x)
print('\nA A^(-1)=')
print(A0*x)
