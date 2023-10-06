#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section  8.6.3                                                     *
%*     filename: ch08pr08.py                                              *
%*     program listing number: 8.8                                        *
%*                                                                        *
%*     This program calculate the determinant of distance matrix for      *
%*     a tree graph using Gaussian elimination with partial pivoting.     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2017.                                    *
%**************************************************************************
"""
import numpy as np

# Set a linear equation
N=10
A=[[0. , 1. , 2. , 3. , 4. , 4. , 3. , 4. , 4. , 5. ],
   [1. , 0. , 1. , 2. , 3. , 3. , 2. , 3. , 3. , 4. ],
   [2. , 1. , 0. , 1. , 2. , 2. , 1. , 2. , 2. , 3. ],
   [3. , 2. , 1. , 0. , 1. , 1. , 2. , 3. , 3. , 4. ],
   [4. , 3. , 2. , 1. , 0. , 2. , 3. , 4. , 4. , 5. ],
   [4. , 3. , 2. , 1. , 2. , 0. , 3. , 4. , 4. , 5. ],
   [3. , 2. , 1. , 2. , 3. , 3. , 0. , 1. , 1. , 2. ],
   [4. , 3. , 2. , 3. , 4. , 4. , 1. , 0. , 2. , 3. ],
   [4. , 3. , 2. , 3. , 4. , 4. , 1. , 2. , 0. , 1. ],
   [5. , 4. , 3. , 4. , 5. , 5. , 2. , 3. , 1. , 0. ]]
A=np.matrix(A)
P=np.matrix(np.identity(N),dtype=int)   # permutation matrix
b=np.matrix(np.zeros(N)).transpose()

S=np.zeros(N)
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
        for i in range(n,N):
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

p=P.sum()
D=(-1)**p
for i in range(0,N):
    D=D*A[i,i]

D_GP=-(N-1)*(-2)**(N-2);
print('Gaussian Elimination: {0:f}'.format(D))
print('       Graham-Pollak: {0:d}'.format(D_GP))



