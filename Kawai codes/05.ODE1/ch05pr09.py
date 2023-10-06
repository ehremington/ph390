#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 5.5.3                                                      *
%*     filename: ch05pr09.m                                               *
%*     program listing number: 5.9                                        *
%*                                                                        *
%*     This program solves the synchronization of two phase oscillators   *
%*     using 2nd-order Runge-Kutta methods.                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

omega1=1.0;  omega2=2.2

# Control parameters
tmax=20; N=2000; h=tmax/N

# Set arrays
t=np.linspace(0,tmax,N+1)
theta1=np.zeros(N+1)
theta2=np.zeros(N+1)

# initial conditions

theta1[0]=np.pi; theta2[0]=0.0

# 2nd-order Runge-Kutta method
for n in range(0,N):
    k1 = omega1 + np.sin(theta2[n]-theta1[n])
    l1 = omega2 - np.sin(theta2[n]-theta1[n])
    mid1 = theta1[n]+k1*h/2.0
    mid2 = theta2[n]+l1*h/2.0
    k2 = omega1 + np.sin(mid2-mid1)
    l2 = omega2 - np.sin(mid2-mid1);
    theta1[n+1]=theta1[n]+k2*h;
    theta2[n+1]=theta2[n]+l2*h;


# plot the trajectories of oscillators
plt.ioff()
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(t,np.sin(theta1),'-b',label='$\theta_1$')
plt.plot(t,np.sin(theta2),'-r',label="$\theta_2$")
plt.xlabel('t')
plt.ylabel('$\sin(\theta)$')
plt.legend(loc=1)

# plot the phase difference
plt.subplot(1,2,2)
plt.plot(t,theta1-theta2,'-k')
plt.xlabel('t');
plt.ylabel('$theta_1-theta_2$')
plt.show()
