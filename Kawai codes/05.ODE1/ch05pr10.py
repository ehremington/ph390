#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
%**************************************************************************
%*     Section 5.5.4                                                      *
%*     filename: ch05pr10.py                                              *
%*     program listing number: 5.10                                       *
%*                                                                        *
%*     This program solves Newton equation for simple harmonic oscillator *
%*     using Runge-Kutta 4th order methods.  Then, it determines the      *
%*     period of the oscillation.                                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/27/2017.                                             *
%**************************************************************************
"""
import numpy as np

# system parameter
omega=1.0; mass=1.0; spring_k=mass*omega**2
E=1.0

# control parameters
h=0.01
Nmax=10

# initial conditions
x0=0.0
v0=np.sqrt(2*(E-spring_k*x0**2/2)/mass)
t0=0.0

N=0

# the first Euler step
x1=x0
x2 = x0 + v0*h - omega**2*x0*h**2/2.0
t=t0
while N<Nmax:
    # Verlet method
    x3=2.0*x2-x1 - omega**2*x2*h**2
    v = (x3-x1)/(2.0*h)
    t=t+h

    # Check if it returned to the tarting point
    if x3-x0>0 and x2-x0< 0 :
        N=N+1
    
    x1=x2
    x2=x3

x=x1-x0
# adjustment of the return time
Fn = -omega**2*x/mass;
if Fn>0 :
    delta = -2.0*x/(v+np.sqrt(v**2-2.0*x*Fn))
else:
    delta = (-v-np.sqrt(v**2-2*x*Fn))/Fn


tau = t+delta
period = tau/N
print('Period: Verlet = {0:7.6f},  Exact = {1:7.6f}'.format(period,2*np.pi))



