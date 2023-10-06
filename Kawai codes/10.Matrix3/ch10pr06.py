#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 10.5.2 - 10.5.4                                            *
%*     filename: ch10pr06.py                                              *
%*     program listing number: 10.6                                       *
%*                                                                        *
%*     This program finds energy and wavefunction of a chain of atoms     *
%*     using the QR algorithm method.                                     *
%*                                                                        *
%*     Uses NUMPY method qr()                                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/18/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

N=10  # The size of the chain molecule
alpha=-2.0
beta=-1.0

# Define the tridiagonal matrix
A=np.zeros((N,N))
for i in range(0,N):
    A[i,i]=alpha
for i in range(0,N-1):
    A[i,i+1]=beta
    A[i+1,i]=beta
A=np.matrix(A)

# Tolerance
tol = 1.0e-8

# Magnitude of total off-diagonal elemens
S=0.0
for i in range(0,N):
    for j in range(0,i):
        S=S+A[i,j]**2

P=np.matrix(np.identity(N))

n=0
while S > tol:
    n+=1
    Q, R=np.linalg.qr(A);  # QR decomposition
    A=R*Q  #Orthogonal transformation
    P=P*Q    #Accumulating the transformation
    # Error evaluation
    S=0.0
    for i in range(0,N):
        for j in range(0,i):
            S=S+A[i,j]**2

print('# of iterations={0:d}'.format(n))
print('Eigenvalues\n')
for n in range(0,N):
    E=-4.*np.sin((N-n)*np.pi/(2*(N+1)))**2
    print('n={0:d} : E={1:f} ({2:f})'.format(n,A[n,n],E))

print('\nTransformation Matrix\n')
print(P)

plt.figure(figsize=(6,5))
plt.plot(np.linspace(0,N-1,N),P[:,0],'-ob',label=r'$\psi_1$')
plt.plot(np.linspace(0,N-1,N),P[:,1],'-og',label=r'$\psi_2$')
plt.plot(np.linspace(0,N-1,N),P[:,2],'-or',label=r'$\psi_3$')

plt.xlabel('x',fontsize=14)
plt.ylabel(r'$\psi$',fontsize=14)
plt.legend(loc=4)
plt.show()