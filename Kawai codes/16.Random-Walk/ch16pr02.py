#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Exercise 16.2                                                      *
%*     filename: ch16pr02.m                                               *
%*     program listing number: 16.2                                       *
%*                                                                        *
%*     This program simulates the one-dimensional persistent random walk. *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/04/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# parameters
p=0.75   # persistent jump probability
M=100000  # number of particles
N=1000    # number of steps

mu=np.zeros(N+1)
sigma2=np.zeros(N+1)
t=np.linspace(0,N,N+1)

# initial states
x=np.zeros(M)  # reset the trajectories
mu[0]=0.
sigma2[0]=0.

d=np.random.choice([-1,1],M)
x=x+d # unbiased initial jump
mu[1]=np.sum(x)/M        # mean
sigma2[1]=np.sum(x**2)/M-mu[1]**2  #variance

for i in range(1,N+1):
    r=np.random.rand(M)
    k=r>p  # direction is reversed.
    d[k]=-d[k]
    x=x+np.sign(d)        # jump
    mu[i]=np.sum(x)/M;       # mean
    sigma2[i]=np.sum(x**2)/M  # ariance

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(t,mu,'-b',label=r'$\mu$',linewidth=2)
plt.plot(t,sigma2,'-r',label=r'$\sigma^2$',linewidth=2)
plt.plot(t,t,'--k',label='$\sigma^2$ (Normal RW)')
plt.legend(loc=2)
plt.xlabel('steps')
plt.ylabel('moments')

plt.subplot(1,2,2)
n, X, Y = plt.hist(x,51,normed=True,label='persistent RW')
dX=X[2]-X[1]
Z=1.0/np.sqrt(2.0*np.pi*N)*np.exp(-X**2/(2.0*N))
plt.plot(X,Z,'-r',label='normal RW')
plt.xlabel('x')
plt.ylabel('probability distribution')
plt.legend(loc=1)
plt.show()
