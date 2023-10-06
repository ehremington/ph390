#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 15.1                                                       *
%*     filename: ch15pr01.py                                              *
%*     program listing number: 16.1                                       *
%*                                                                        *
%*     This program simulate a dice using a psueo random number           *
%*     generator.                                                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# parapmeters for random number generators
a=np.int64(16807)
b=np.int64(0) 
c=np.int64(2147483647)

# get a seed
x=np.int64(input('Seed='))

# generate uniform random numbers
N=6000
r=np.zeros(N)
for i in range(0,N):
    x=np.mod(a*x,c)
    r[i]=np.float(x)/np.float(c)
 
# statistics of virtual die
P=np.array([0,0,0,0,0,0])
D=np.array([1,2,3,4,5,6])
for i in range(0,N):
    n=np.int(np.ceil(6.0*r[i]))-1
    P[n]=P[n]+1;

plt.figure(figsize=(6,5))
plt.bar(D,P,0.9)
plt.plot([0.5,6.5],[N/6,N/6],'--r')
plt.xlabel('face of die',fontsize=14)
plt.ylabel('number of realization',fontsize=14)
plt.show()
