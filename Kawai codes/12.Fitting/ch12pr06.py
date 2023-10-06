#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 12.6                                                       *
%*     filename: ch12pr06.py                                              *
%*     program listing number: 12.6                                       *
%*                                                                        *
%*     This program interpolates 11-point data with the quadratic         *
%*     regression method.                                                 *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# Generate a noisy data set
N=13
x=np.zeros(N)
y=np.zeros(N)
s=np.zeros(N)
sm=5.0
for i in range(0,N):
    x[i]=i-6.+np.random.uniform(-0.2,0.2)
    s[i]=np.random.normal(0.0,sm/2.)+sm
    y[i]=-2*x[i]**2+s[i]*np.random.uniform(0.2,0.9)

M=3; # number of parameters
J=np.matrix(np.zeros((N,M)))
b=np.matrix(np.zeros(N)).transpose()

# Construct Jacobian matrix and a vector
for i in range(0,N):
    for j in range(0,M):
        J[i,j]=x[i]**j/s[i]

    b[i]=y[i]/s[i]

# Solve the linear equation
A=J.transpose()*J
c=J.transpose()*b
lam=np.linalg.solve(A,c)

# constructe the fitted curve
K=121
z=np.linspace(-6.0,6.0,K)
f=lam[0,0]+lam[1,0]*z+lam[2,0]*z**2

plt.figure(figsize=(6,5))
plt.errorbar(x,y,yerr=s,fmt='ok',label='Data')
plt.plot(z,f,'-r',label='Fit')
plt.xlabel('x',fontsize=14)
plt.ylabel('f(x)',fontsize=14)
plt.legend(loc=1)
plt.show()

