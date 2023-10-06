#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 10.5.1                                                     *
%*     filename: ch10pr05.py                                              *
%*     program listing number: 10.5                                       *
%*                                                                        *
%*     This program finds eigenmodes of cpoupled harmonic oscillators.    *
%*     The higest and lowest frequencies are obgtained by the power       *
%*     method and the remaining frequency is computed by the inverse      *
%*     iteration method.                                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/18/2017.                                    *
%**************************************************************************
"""
import numpy as np

# system parameters
k1=2.0; k2=4.0; k3=4.0; k4=2.0
m1=2.0; m2=4.0; m3=3.0
K=np.matrix([[k1+k2, -k2, 0.0],[-k2,k2+k3,-k3],[0.0,-k3,k3+k4]])
Minv=np.matrix([[1.0/m1,0.0,0.0],[0.0,1.0/m2,0.0],[0.0,0.0,1/m3]])
A0=Minv*K
A=np.matrix(np.zeros((3,3)))
u=np.matrix(np.zeros((3,3)))
eig=np.zeros(3)

tol = 1e-7 # tolerance

# Find the largest/smallest eigenvalues
# by the power method.
for i in (0,1,2):

    if i==0:
        A[:,:]=A0[:,:]  # for largest eigenvalue
    else:
        A=np.linalg.inv(A0) # for smallest

    found=False
    x=np.matrix(np.random.rand(3)).transpose()    # initial guess
    u0=x/np.linalg.norm(x)  # normalization

# power method iteration
    n=1
    while not(found):
        x=A*x   # update x
        u1=x/np.linalg.norm(x)   # normalization
        err=np.linalg.norm(u1-u0)   # error
        if err < tol:
            found=True
        else:
            u0=u1
            n+=1

        if i==0:
            eig[0]=u1.transpose()*A*u1 # largest eigenvalue
            u[:,0]=u1
        else:
            eig[2]=1.0/(u1.transpose()*A*u1) # smallest
            u[:,2]=u1

# The other eigenvalue by the inverse 
# iterative method.
# guess=middle between the largest and
# the lowest.

u0=np.matrix(np.random.rand(3)).transpose() # random vector
u0=u0/np.linalg.norm(u0)
q = (eig[0]+eig[2])/2.0 
A=A0-q*np.identity(3)  # shifted matrix

found=False
n=1

# inverse iteration method
while not(found):
    x=np.linalg.solve(A,u0)
    u1=x/np.linalg.norm(x)
    if u1[0]<0:  # correct the phase.
        u1=-u1

    err=np.linalg.norm(u1-u0) 
    if err < tol:
        found=True
    else:
        u0=u1
        n+=1

# eigenvalue
eig[1]=u1.transpose()*A0*u1
u[:,1]=u1

# output
for i in range(0,3):
    print('\nFrequency={0:f}'.format(np.sqrt(eig[i])))
    print('Eigenvector')
    print(u[:,i])

