#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 12.1                                                       *
%*     filename: ch12pr01.py                                              *
%*     program listing number: 12.1                                       *
%*                                                                        *
%*     This program interpolates 11-point data with linear spline.        *
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
N=F.size # Num er of data points
h=np.zeros(N-1)
for j in range(0,N-1):
    h[j]=X[j+1]-X[j]

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
        y[n]=(1.0-t)*F[i]+t*F[i+1]
        n+=1

x[n]=X[N-1]
y[n]=F[N-1]
z=np.sin(x)*np.exp(-0.2*x)

n+=1
plt.figure(figsize=(6,5))
plt.plot(x[0:n],y[0:n],'-r',label='Spline')
plt.plot(X,F,'ob',label='Data')
plt.plot(x[0:n],z[0:n],'--k',label='Source')
plt.xlabel('x',fontsize=14)
plt.ylabel('f(x)',fontsize=14)
plt.legend(loc=1)
plt.show()
