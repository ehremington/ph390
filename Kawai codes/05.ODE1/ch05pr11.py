#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 5.5.5                                                      *
%*     filename: ch05pr11.py                                              *
%*     program listing number: 5.11                                       *
%*                                                                        *
%*     This program finds the trajectory of a pendulum using              *
%*     Euler and Verlet methods.    Euler method shows its numerical      *
%*    instability and the trajectory diverges.                            *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/27/2017.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# system parameters
mass=1.0; L=1.0; g=9.8; Omega=np.sqrt(g/L); I=mass*L**2

# control parameters
tmax=50; N=5000; h=tmax/np.float(N)

theta1=np.zeros(N+1)
theta2=np.zeros(N+2)
omega1=np.zeros(N+1)
omega2=np.zeros(N+1)
E1=np.zeros(N+1)
E2=np.zeros(N+1)
t=np.linspace(0.0,tmax,N+1)

# initial conditions
theta1[0]=0.5; omega1[0]=0.0
theta2[0]=0.5; omega2[0]=0.0
E1[0] = - mass*g*L*np.cos(theta1[0])
E2[0] = - mass*g*L*np.cos(theta2[0])


# Euler method
for i in range(0,N):
    omega1[i+1] = omega1[i] - Omega**2 * np.sin(theta1[i]) * h
    theta1[i+1] = theta1[i] + omega1[i]*h
    E1[i+1] = I/2.0 * omega1[i+1]**2 - mass*g*L*np.cos(theta1[i+1])


# Verlet method
theta2[1] = theta2[0] + omega2[0]*h- Omega**2 * np.sin(theta2[0])*h**2/2.0
for i in range(1,N+1):
    theta2[i+1]=2*theta2[i]-theta2[i-1]-Omega**2 * np.sin(theta2[i])*h**2.0
    omega2[i] = (theta2[i+1]-theta2[i-1])/(2.0*h)
    E2[i]=I/2 * omega2[i]**2 - mass*g*L*np.cos(theta2[i])

plt.ioff()
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(t,theta1,'-b',label='Euler')
plt.plot(t,theta2[0:N+1],'-r',label='Verlet')
plt.xlabel('t')
plt.ylabel('$\theta$')
plt.legend(loc=3)

plt.subplot(1,2,2)
plt.plot(t,E1,'-b',label='Euler')
plt.plot(t,E2,'-r',label='Verlet')
plt.xlabel('t')
plt.ylabel('Energy')
plt.legend(loc=2)
plt.show()