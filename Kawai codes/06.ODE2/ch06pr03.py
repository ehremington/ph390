#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section  6.3.1                                                     *
%*     filename: ch06pr03.m                                               *
%*     program listing number: 6.3                                        *
%*                                                                        *
%*     This program finds the wave function of freely falling particle.   *
%*     to reach height yf in travel time tf.                              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import airy 

# define w(x) in Numerov method
def W(x):
    return -x

# control parameters
xmax = 5.0
xmin = -15.0

# Integrating from xmax to 0
N = 200
h = (xmin-xmax)/np.float(N)
x = np.linspace(xmax,xmin,N+1)
y = np.zeros(N+1)
w = np.zeros(N+1)

# initial conditions
y[0] = 0
w[0] = W(x[0])

# we guess next value
y[1] = y[0]+0.1
w[1] = W(x[1]);

# shoot left by the Numerov method
for n in range(1,N):
  w[n+1] = W(x[n+1]);
  y[n+1] = 2.0*(1.0-5.0*h**2*w[n]/12.0)*y[n] - (1.0+h**2*w[n-1]/12.0)*y[n-1]
  y[n+1] = y[n+1]/(1.0+h**2*w[n+1]/12.0)

# normalization
N0 = np.int(-xmax/h) # find the location of x=0
y[:] = y[:]/y[N0]


plt.figure(figsize=(6,5))
plt.plot(x,y,'-r',label="Numerov",linewidth=2.5)
plt.plot(x,airy(x)[0]/airy(0)[0],'-b',label="Exact")
plt.xlabel('x')
plt.ylabel(r'$\psi(x)$')
plt.legend(loc=4)
plt.plot([-15, 5],[0,0],'--k',[0,0],[-1.5,2],'--k')
plt.show()
