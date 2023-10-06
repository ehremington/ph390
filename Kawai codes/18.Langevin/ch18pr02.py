#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Exercise 18.2                                                      *
%*     filename: ch18pr02.py                                              *
%*     program listing number: 18.2                                       *
%*                                                                        *
%*     This program simulates the two-dimensional free continuous-time    *
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
gamma=0.1  # frictional constant
T=1.0      # temperature  
D=T/gamma  # strength of noise

# control parameter
tau=100# total time
dt=0.005 # time step
ds=np.sqrt(2.*D*dt) # Wiener step size
N=np.int(tau/dt)  # number of steps
t = np.linspace(0.0,N*dt,N+1)  # time grid
M=1000 # number of samplings


# reset the counters
sumx1=np.zeros(N+1)
sumx2=np.zeros(N+1)
sumy1=np.zeros(N+1)
sumy2=np.zeros(N+1)
x=np.zeros(N+1)
y=np.zeros(N+1)


# sampling loop begins
for j in range(0,M):
    
    x[0]=0; y[0]=0  # initial position
    
    # Normally ditributed random numbers by the Box-Muller method
    g=np.random.randn(N,2)

    # integration of the Langevin equation
    for i in range(0,N):
        x[i+1] = x[i]+g[i,0]*ds
        y[i+1] = y[i]+g[i,1]*ds
    
    # record the data for statistical analysis
    sumx1 = sumx1+x
    sumy1 = sumy1+y
    sumx2 = sumx2+x**2
    sumy2 = sumy2+y**2

# plot a sample trajectory
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(x,y,'-b')
plt.plot(x[0],y[0],'.r')
plt.axis('equal')
plt.xlabel('x',fontsize=14)
plt.ylabel('y',fontsize=14)
plt.show()


# statistical analysis
plt.subplot(1,2,2)
xmean=sumx1/M   # <x>
ymean=sumy1/M   # <y>
sumx2=sumx2/M   # <x^2>
sumy2=sumy2/M   # <y^2>
xvar=sumx2-xmean**2  # variance_x
yvar=sumy2-ymean**2  # variance_y
theory=2.0*T/gamma*t  # theoretical variance
plt.plot(t,xmean,'-b',label=r'$\langle x \rangle')
plt.plot(t,ymean,'-b',label=r'$\langle y \rangle')
plt.plot(t,xvar,'-r',label=r'$\sigma_x^2$')
plt.plot(t,yvar,'-r',label=r'$\sigma_y^2$')
plt.plot(t,theory,'--r',label=r'$\sigma^2_{ex}$')
plt.xlabel('t',fontsize=14)
plt,ylabel('moments',fontsize=14)
plt.legend(loc=2)
plt.show()
