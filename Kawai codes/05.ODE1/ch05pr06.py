#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 5.6                                                        *
%*     filename: ch05pr06.py                                              *
%*     program listing number: 5.6                                        *
%*                                                                        *
%*     This program solves Newton equation for simple harmonic oscillator *
%*     using Verlet method.                                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/22/2017.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# system parameter
omega=1;

# control parameters
tmax=8.0*np.pi/omega
N=500
h=tmax/N

t=np.linspace(0,tmax,N+1)
x=np.zeros(N+1)
v=np.zeros(N+1)

# exact soution
x_ex=np.cos(omega*t)

# initial conditions
x[0]=1.0
v[0]=0.0

# the first Euler step
x[1] = x[0] + v[0]*h - omega**2*x[0]*h**2/2.0

for n in range(1,N):
    # Verlet method
    x[n+1]=2*x[n]-x[n-1] - omega**2*x[n]*h**2;


# plot trajectories
plt.ioff()
plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.plot(t,x,'ob',label='Verlet')
plt.plot(t,x_ex,'-r',label='Exact')
plt.xlabel('t')
plt.ylabel('displacement')
plt.legend(loc=1)

# plot absolute error
plt.subplot(1,2,2)
plt.semilogy(t,abs(x-x_ex),'-k')
plt.xlabel('t');
plt.ylabel('absolute error')
plt.show()
