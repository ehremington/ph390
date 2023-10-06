#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 5.4                                                        *
%*     filename: ch05pr04.py                                              *
%*     program listing number: 5.4                                        *
%*                                                                        *
%*     This program solves Newton equation for interacting two cars       *
%*     using Runge-Kutta 2nd order methods.                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/22/2017.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# Control parameters
tmax=8; N=400; h=tmax/N;

t=np.linspace(0,tmax,N+1)
v1=np.zeros(N+1)
v2=np.zeros(N+1)

#initial conditions
v1[0]=1.2; v2[0]=1.0

# 2nd-order Runge-Kutta method
for n in range(0,N):
    k1 =   v2[n]-v1[n]
    l1 = -(v2[n]-v1[n])
    mid1 = v1[n]+k1*h/2
    mid2 = v2[n]+l1*h/2
    k2 =   mid2-mid1
    l2 = -(mid2-mid1)
    v1[n+1]=v1[n]+k2*h
    v2[n+1]=v2[n]+l2*h

plt.ioff()
plt.figure(figsize=(12,5))

plt.subplot(1,2,1) 
plt.plot(t,v1,'-b',label="$v_1$")
plt.plot(t,v2,'-r',label="$v_2$")
plt.xlabel('t')
plt.ylabel('velocity')
plt.legend(loc=1)

plt.subplot(1,2,2)
plt.plot(t,v1-v2,'-k')
plt.xlabel('t')
plt.ylabel("$v_1-v_2$")
plt.show()

