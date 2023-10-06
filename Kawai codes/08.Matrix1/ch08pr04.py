#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  8.4                                                       *
%*     filename: ch08pr04.pu                                              *
%*     program listing number: 8.4                                        *
%*                                                                        *
%*     This program solves a simple linear equation with the Gaussian     *
%*     elimination with poartial pivoting and backsubstitution methods.   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2017.                                    *
%**************************************************************************
"""

# Set a linear equation
N=3;
A=np.matrix([[3.,-1.,4.],[2.,0.,-1.],[0.,3.,2.]])
b=np.matrix.transpose(np.matrix([2.,-1.,3.]))
x=np.matrix.transpose(np.matrix(np.zeros(N)))
# permutation matrix must be initially an identity matrix
P=np.matrix(np.identity(3,dtype=int))  
S=np.zeros(3)

# Find scale factors
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
    if j != n :
        for i in range(n,N) :
            TMP=A[n,i]
            A[n,i]=A[j,i]
            A[j,i]=TMP

        TMP=np.asscalar(b[n])
        b[n]=b[j]
        b[j]=TMP
        # Record the permutation
        P[n,n]=P[j,j]=0
        P[n,j]=P[j,n]=1

    # Gaussian elimination
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
print('\nP=\n')
print(P)
print('\nx=\n')
print(x)
