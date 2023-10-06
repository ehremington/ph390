#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Exercise 18.1                                                      *
%*     filename: ch18pr01.py                                              *
%*     program listing number: 18.1                                       *
%*                                                                        *
%*     This program simulates the one-dimensional free continuous-time    *
%*     random walk (Wiener process) using the Heun method.                *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/25/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# system parameters
gamma=0.1  # friction coefficient
T=1.0      # temperature

# control parameter
tau=100  # duration
dt=0.005 # time step
D=T/gamma    # diffusion constant
ds=np.sqrt(2.*D*dt)  # time step for Wiener process
N=np.int(tau/dt)  # number of steps
M=1000  # number of samplings
t = np.linspace(0.0,N*dt,N+1)  # time grid

# create arrays
sumx1=np.zeros(N+1)
sumx2=np.zeros(N+1)
x=np.zeros(N+1)
plt.figure()
theory=2.0*T/gamma*t    # theoretical variance
plt.plot([t[0],t[N]],[0,0],'--r')
plt.plot(t,sqrt(theory),'-r')
plt.plot(t,-sqrt(theory),'-r')

# loop over M samplings
for j in range(0,M):
    # initial condition
    x[0]=0

    # random force with Gaussian stribution
    g=np.random.randn(N)  

    # solving Langevin equation
    for i in range(0,N):
        x[i+1] = x[i]+g[i]*ds
    if j<10: 
        plt.plot(t,x,)
        plt.pause(0.0001)
        
    # save the data for statistical analysis
    sumx1 = sumx1+x
    sumx2 = sumx2+x**2

plt.xlabel('t',fontsize=14)
plt.ylabel('x(t)',fontsize=14)
plt.show()

# statistical analysis
xmean=sumx1/M    # <x>
sumx2=sumx2/M    # <x^2>
xvar=sumx2-xmean**2   # variance


plt.figure()
plt.plot(t,xmean,'-b',label='mean')
plt.plot(t,xvar,'-g',label='variance')
plt.plot(t,theory,'--g',label='var_exact')
plt.xlabel('t',fontsize=14)
plt.ylabel('mean, vriance',fontsize=14)
plt.legend(loc=2)
plt.show()
