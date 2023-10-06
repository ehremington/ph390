#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 10.5                                                       *
%*     filename: ch10pr03.m                                               *
%*     program listing number: 10.3                                       *
%*                                                                        *
%*     This program finds eigenvalues and eigenvectos of a symmetric      *
%*     matrix using the Jacobi transformation method.                     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2015.                                    *
%**************************************************************************
"""
import numpy as np

# Define the matrix
A=np.matrix([[1.,-4.,2.],[-4.,1.,-2.],[2.,-2.,-2.]])

# Tolerance
tol=1.0e-4

# Evaluate error
S=0.0
for i in range(0,3):
    for j in range(i+1,3):
        S=S+abs(A[i,j])

# Initial transformation matrix
P=np.matrix(np.identity(3))

while S > tol:
    for i in range(0,3):
        for j in range(i+1,3):
            if A[i,j] != 0:
                # Jacobian rotation
                if A[j,j]==A[i,i]:
                    c=-1./np.sqrt(2.)
                    s=-1./np.sqrt(2.)
                else:
                    beta=(A[j,j]-A[i,i])/(2.0*A[i,j])
                    t=np.sign(beta)/(np.abs(beta)+np.sqrt(beta**2+1.0))
                    c=1./np.sqrt(t**2+1.0)
                    s= t/np.sqrt(t**2+1.0)

                r=s/(1.0+c)
                ai = c**2*A[i,i]+s**2*A[j,j]-2.0*s*c*A[i,j]
                aj = s**2*A[i,i]+c**2*A[j,j]+2.0*s*c*A[i,j]
                A[i,i]=ai
                A[j,j]=aj
                # Transformation matrix
                Q=np.matrix(np.identity(3))
                Q[i,i]=c
                Q[j,j]=c
                Q[i,j]=s
                Q[j,i]=-s
                P=P*Q

                for k in range(0,3):
                    if k!=i and k!=j:
                        aki=c*A[k,i]-s*A[k,j]
                        akj=c*A[k,j]+s*A[k,i]
                        A[k,i]=aki
                        A[i,k]=aki
                        A[k,j]=akj
                        A[j,k]=akj

                A[i,j]=0.0
                A[j,i]=0.0

    # Evaluate error
    S=0.0
    for i in range(0,3):
        for j in range(i+1,3):
            S=S+abs(A[i,j])

print('\nTransformed Matrix\n')
print(A)

print('\nTransformation Matrix\n')
print(P)
