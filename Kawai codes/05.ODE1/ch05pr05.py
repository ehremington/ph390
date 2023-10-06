#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 5.5                                                        *
%*     filename: ch05pr05.py                                              *
%*     program listing number: 5.5                                        *
%*                                                                        *
%*     This program solves Newton equation for simple harmonic oscillator *
%*     using Runge-Kutta 4th order methods.                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/22/2017.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# system parameter
omega=1.0

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

for n in range(0,N):
    # 4th-order Runge-Kutta
    kv1=-omega**2*x[n];
    kx1=v[n];
    
    v_mid = v[n]+kv1*h/2.0
    x_mid = x[n]+kx1*h/2.0
    kv2 = -omega**2*x_mid
    kx2 = v_mid
      
    v_mid = v[n]+kv2*h/2.0
    x_mid = x[n]+kx2*h/2.0
    kv3 = -omega**2*x_mid
    kx3 = v_mid

    v_end = v[n]+kv3*h
    x_end = x[n]+kx3*h
    kv4 = -omega**2*x_end
    kx4 = v_end
    
    v[n+1]=v[n]+(kv1+2.0*(kv2+kv3)+kv4)*h/6.0
    x[n+1]=x[n]+(kx1+2.0*(kx2+kx3)+kx4)*h/6.0

# plot trajectories
plt.ioff()
plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.plot(t,x,'ob',label='RK4')
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

