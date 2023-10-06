#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 5.5.1                                                      *
%*     filename: ch05pr07.py                                              *
%*     program listing number: 5.7                                        *
%*                                                                        *
%*     This program solves the Brusselator model                          *
%*     using Runge-Kutta 4th order methods.                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/22/2017.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# fixed parameter
a=1.0

# control parameter
b=np.float(input('Enter value for b [1.5-2.5] ='))

# duration 
tmax=100
# number of integration steps
N=2000
# step size
h=tmax/N

t=np.linspace(0,tmax,N+1)
x=np.zeros(N+1)
y=np.zeros(N+1)

#initial conditions
x[0]=1.0
y[0]=1.0

for n in range(0,N):
    # 4th-order Runge-Kutta
    kx1=a-(b+1)*x[n]+x[n]**2*y[n]
    ky1=b*x[n]-x[n]**2*y[n]
    
    x_mid = x[n]+kx1*h/2.0
    y_mid = y[n]+ky1*h/2.0
    kx2 = a - (b+1.0)*x_mid+x_mid**2*y_mid
    ky2 = b*x_mid-x_mid**2*y_mid
      
    x_mid = x[n]+kx2*h/2.0
    y_mid = y[n]+ky2*h/2.0
    kx3 = a - (b+1.0)*x_mid+x_mid**2*y_mid
    ky3 = b*x_mid-x_mid**2*y_mid

    x_end = x[n]+kx3*h
    y_end = y[n]+ky3*h
    kx4 = a - (b+1.0)*x_end+x_end**2*y_end
    ky4 = b*x_end-x_end**2*y_end
    
    x[n+1]=x[n]+(kx1+2.0*(kx2+kx3)+kx4)*h/6.0
    y[n+1]=y[n]+(ky1+2.0*(ky2+ky3)+ky4)*h/6.0

# plot individual trajectories
plt.ioff()
plt.figure(figsize=(12,5))

plt.subplot(1,2,1);
plt.plot(t,x,'-b',label='x')
plt.plot(t,y,'-r',label='y')
plt.xlabel('t');
plt.ylabel('concentration');
plt.legend(loc=1)

# plot phase trajectory
plt.subplot(1,2,2)
plt.plot(x,y,'-b')
plt.xlabel('x')
plt.ylabel('y');
plt.show()
