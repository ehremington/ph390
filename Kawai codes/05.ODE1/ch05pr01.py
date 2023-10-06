#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 5.1                                                        *
%*     filename: ch05pr01.py                                              *
%*     program listing number: 5.1                                        *
%*                                                                        *
%*     This program solves Newton equation for a falling object           *
%*     using Euler and predictor-corrector methods.                       *
%*        m = mass of the object                                          *
%*        g = acceleration due to gravity                                 *
%*        gamma = frictional coefficient                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/20/2017.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# system parameters
gamma=1.0
g=9.8
m=1.0

# integration parameters
tmax=10  # maximum time
N=1000   # maximum steps
h=tmax/N # time step

# set arrays
v_ex=np.zeros(N+1)
v_eu=np.zeros(N+1)
v_pc=np.zeros(N+1)
t=np.linspace(0,N,N+1)*h

# intial condition
v_ex[0]=0.0
v_eu[0]=0.0
v_pc[0]=0.0

for i in range(0,N):
    # Euler method
    F_eu = -gamma*v_eu[i]/m - g 
    v_eu[i+1] = v_eu[i] + F_eu*h
    
    # Predictor-Corrector method
    F_pc = -gamma*v_pc[i]/m - g
    v_pc[i+1] = v_pc[i] + F_pc*h;  # predictor
    F_pc = -gamma/m*(v_pc[i]+v_pc[i+1])/2 - g
    v_pc[i+1] = v_pc[i] + F_pc*h   # corrector
    
    # Exact solution
    v_ex[i+1] = m*g/gamma*(np.exp(-gamma*t[i+1])-1)

plt.ioff()
plt.figure(figsize=(12,5))

# Plot the solutions
plt.subplot(1,2,1);
plt.plot(t,v_eu,'-b',label='Euler')
plt.plot(t,v_pc,'-r',label='Predictor-Corrector')
plt.plot(t,v_ex,'-k',label='Exact')
plt.xlabel('t')
plt.ylabel('v(t)')
plt.legend(loc=1)

# Plot the absolute errors
plt.subplot(1,2,2)
plt.semilogy(t,abs(v_eu-v_ex),'-b',label='Euler')
plt.semilogy(t,abs(v_pc-v_ex),'-r',label='Predictor-Corrector')

plt.xlabel('t')
plt.ylabel('absolute error')
plt.legend(loc=3)
plt.show()