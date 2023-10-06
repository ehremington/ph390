#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 15.2                                                       *
%*     filename: ch15pr02.m                                               *
%*     program listing number: 15.2                                       *
%*                                                                        *
%*     This program  evaluate the value of pi using the Monte Carlo       *
%*     integeration of a circle.                                          *
%*     Uses: numpy random package                                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# random points on a square
N=100000
x=np.random.uniform(-1.0,1.0,N)
y=np.random.uniform(-1.0,1.0,N)

# points inside the circle
hit=0.0
M=np.int(N/100)
PI=np.zeros(M)
i=0
for n in range(0,N):
    if x[n]**2+y[n]**2 < 1.0:
        hit+=1.0
        
    if n>0 and np.mod(n,100)==0: # evaluate at every 100 
        PI[i]=hit/n*4.0 # estimate of pi
        i+=1

plt.figure(figsize=(6,5))
T=np.linspace(100,N,M)
plt.plot(T,PI)
plt.plot([0, N], [np.pi,np.pi],'--r')
plt.xlabel('# of sampling',fontsize=14)
plt.ylabel(r'$\pi$ by Monte Carlo integration',fontsize=14)
plt.axis([0, N, np.pi*0.9, np.pi*1.1])
plt.show()