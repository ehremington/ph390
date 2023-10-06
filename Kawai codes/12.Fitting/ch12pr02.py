#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 12.2                                                       *
%*     filename: ch12pr02.py                                              *
%*     program listing number: 12.2                                       *
%*                                                                        *
%*     This program interpolates 11-point data with cubic spline.         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# data to be fitted
F=[0.0000, 0.6889, 0.6095, 0.0774, -0.3401, -0.3528,\
   -0.0842, 0.1620, 0.1997, 0.0681, -0.0736]
X=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.]
F=np.array(F)
X=np.array(X)
N=F.size
h=np.zeros(N-1)
for j in range(0,N-1):
    h[j]=X[j+1]-X[j]

G=np.zeros(N-2)
for j in range(0,N-2):
    G[j]=3*((F[j+2]-F[j+1])/h[j+1]-(F[j+1]-F[j])/h[j])

A=np.zeros((N-2,N-2))
for j in range(0,N-2):
    A[j,j]=(h[j+1]+h[j])/2.0

for j in range(0,N-3):
    A[j,j+1]=h[j+1]
    A[j+1,j]=h[j+1]

P=np.linalg.solve(A,G)

Q=np.zeros(N)
for j in range(0,N-2):
    Q[j+1]=P[j]

Q[0]=0.0
Q[N-1]=0.0

M=10  # Number of interpolation points between data points
dt=1.0/M
T=np.linspace(0.0,dt*(M-1),M) # linear interpolation
x=np.zeros(N*M)
y=np.zeros(N*M)
n=0
for i in range(0,N-1):
    # linear interpolation between two adjacent data points
    for t in T:
        x[n]=t*h[i]+X[i]
        y[n]=h[i]**2/6.0 * Q[i]   * t*(t+1.0)*(t-1.0) \
            -h[i]**2/6.0 * Q[i+1] * t*(t-1.0)*(t-2.0) \
            +F[i+1]*t + (1.0-t)*F[i]
        n+=1

z=np.sin(x[0:n])*np.exp(-0.2*x[0:n])

plt.figure(figsize=(6,5))
plt.plot(x[0:n],y[0:n],'-r',label='Spline')
plt.plot(X,F,'ob',label='Data')
plt.plot(x[0:n],z[0:n],'--k',label='Source')
plt.xlabel('x',fontsize=14)
plt.ylabel('f(x)',fontsize=14)
plt.legend(loc=1)
plt.show()
