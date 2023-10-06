#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Exercise 17.1                                                      *
%*     filename: ch17pr01.m                                               *
%*     program listing number: 17.1                                       *
%*                                                                        *
%*     This program finds the velocity ditribution of one-dimensional     *
%*     ideal gas using Metropolis algorithm.                              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/11/2017.                                    *
%**************************************************************************
"""

import numpy as np
import matplotlib.pyplot as plt

N0=10000 # number of thermalization steps
N=200000; # number of samples
kT=3.0 # Temperature times Boltzmann constant
m=5.0 # mass of particle

v=np.zeros(N+N0)
dv=0.1    # maximum jump in velocity
v0=np.sqrt(kT/m)   # thermal speed
# initial velocity (uniform random beteen -v0 and +v0)
v[0]=2.*v0*(np.random.rand(1)-0.5)  

for i in range(0,N+N0-1):
    found = False
    while not(found):
        u = v[i] + dv*(2.0*np.random.rand(1)-1.0)   # candidate
        dE = m/2.0*(u**2-v[i]**2)           # energy change
        if np.exp(-dE/kT)>np.random.rand(1):    # Metropolis condition
            v[i+1]=u                         # accept change ()
            found = True

# theoretical distribution (Maxwell)
K=41
plt.close()
plt.figure(figsize=(6,5))
n, bins, patches = plt.hist(v[N0:N0+N],K,normed=1,label='Monte Carlo')
w=np.linspace(bins[0],bins[-1],101)
y=1.0/np.sqrt(2.*np.pi*kT/m) * np.exp(-m*w**2/(2.*kT))
plt.plot(w,y,'-r',label='Maxwell')
plt.legend(loc=1)
plt.xlabel('v')
plt.ylabel('p(v)')
plt.show()


