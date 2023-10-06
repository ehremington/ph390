#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 12.3                                                       *
%*     filename: ch12pr03.py                                              *
%*     program listing number: 12.3                                       *
%*                                                                        *
%*     This program interpolates 11-point data with the Vandermonde       *
%*     matrix.                                                            *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/24/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

F=[0.0000, 0.6889, 0.6095, 0.0774, -0.3401, -0.3528,\
   -0.0842, 0.1620, 0.1997, 0.0681, -0.0736]
X=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.]
F=np.array(F)
X=np.array(X)
N=F.size

# construction of the Vandermonde matrix
x = np.zeros((N,N))
x[:,0]=1.0
for n in range(1,N):
    x[:,n]=X[:]**n

# solve the linear equation
# using Gaussian elimination
a=np.linalg.solve(x,F)

# evaluate the function value
# between the sampling points.

M=101
z=np.linspace(0.0,X[N-1],M)
y=np.zeros(M)

for j in range(0,M):
    y[j]=a[0]
    for i in range(1,N):
        y[j]=y[j]+a[i]*z[j]**i

v = np.sin(z)*np.exp(-0.2*z)

plt.figure(figsize=(6,5))
plt.plot(X,F,'ob',label="Raw data")
plt.plot(z,y,'-r',label="Vandermonde")
plt.plot(z,v,'--k',label="Source")
plt.xlabel('x',fontsize=14)
plt.ylabel('f(x)',fontsize=14)
plt.legend(loc=1)
plt.show()
