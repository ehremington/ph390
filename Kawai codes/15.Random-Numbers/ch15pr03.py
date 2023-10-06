#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 16.3                                                       *
%*     filename: ch16pr03.m                                               *
%*     program listing number: 16.3                                       *
%*                                                                        *
%*     This program generates Gaussian distributed random numbers         *
%*     using the Box-Muller method.                                       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/15/2015.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

N=10000000
x=np.random.rand(N,2)
y=np.zeros((N,2))
y[:,0] = np.sqrt(-2.0*np.log(x[:,0]))*np.cos(2.0*np.pi*x[:,1])
y[:,1] = np.sqrt(-2.0*np.log(x[:,0]))*np.sin(2.0*np.pi*x[:,1])
y=y.reshape(2*N)
ymin=min(y)
ymax=max(y)
print('the largest deviation={0:f}, {1:f}\n'.format(ymin,ymax))

K=201
plt.figure(figsize=(6,5))
plt.hist(y,K,normed=1,label='Box-Muller')
xmin=np.floor(ymin)
xmax=np.ceil(ymax)
X=np.linspace(xmin,xmax,K)
F=np.exp(-X**2/2.)/np.sqrt(2.*np.pi)
plt.plot(X,F,'-r',label='Exact')
plt.xlabel('x',fontsize=14)
plt.ylabel('p(x)',fontsize=14)
plt.legend(loc=1)
plt.show()
