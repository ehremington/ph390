#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Exercise 16.1                                                      *
%*     filename: ch16pr01.py                                              *
%*     program listing number: 16.1                                       *
%*                                                                        *
%*     This program simulates the one-dimensional descrete random walk.   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/15/2014.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

N=1000    # max time (number of steps)
M=100000  # number of particles

x=np.zeros((M,N+1))  # reset the trajectories

# trajectory calculation
for i in range(0,N):
    x[:,i+1]=x[:,i]+np.random.choice([-1,1],M) # random step

# stattistics
mu=np.sum(x,axis=0)/M        # mean
sigma=np.sqrt(np.sum(x**2,axis=0)/M)-mu**2  #variance
t=np.linspace(0,N,N+1)

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(t,mu,'-r',label=r'$\mu$')
plt.plot(t, sigma,'--k',label=r'$+\sigma$')
plt.plot(t,-sigma,'--k',label=r'$-\sigma$')
plt.plot(t,x[0,:],'-b')
plt.plot(t,x[1,:],'-g')
plt.plot(t,x[2,:],'-c')
plt.axis([0, N, -40, 40])
plt.xlabel('steps',fontsize=14)
plt.ylabel('x',fontsize=14)
plt.legend(loc=3)

plt.subplot(1,2,2)
rho=x[:,N]
n, X, Y = plt.hist(rho,51,normed=True,label='simulation')
Y=1.0/np.sqrt(2*np.pi*N)*np.exp(-X**2/(2.*N))
dX=X[2]-X[1]
plt.plot(X,Y,'-r',label='theory')
plt.xlabel('x')
plt.ylabel(r'$\rho(x)$')
plt.legend(loc=1)
plt.show()