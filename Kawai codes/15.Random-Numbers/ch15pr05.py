#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 15.5.1                                                     *
%*     filename: ch15pr05.m                                               *
%*     program listing number: 15.5                                       *
%*                                                                        *
%*     This program evaluate the mean speed of the gas particles          *
%*     in a thermal equilibrium.                                          *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
"""
import numpy as np

# parameters
T=300.             # Temperature in K
k=1.380658e-23    # Boltzman constant in J/K
m=2*1.672623e-27  # H2 mass in kg

# velocity at Maxwell distribution
N=100000000
s=np.sqrt(k*T/m)
v=np.random.normal(0.0,s,[N,3]) # 3 components (vx, vy, vz)

# speed
speed=np.sqrt(v[:,0]**2+v[:,1]**2+v[:,2]**2)
# mean
mean=sum(speed)/N
#theory
exact=2*s*np.sqrt(2/np.pi)
# error
error=np.abs(mean-exact)/exact

print('mean speed = {0:10.5e} (exact={1:10.5e})'.format(mean,exact))
print('relative error = {0:10.5}\n'.format(error))

