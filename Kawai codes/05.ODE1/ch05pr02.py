#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 5.2                                                        *
%*     filename: ch05pr02.m                                               *
%*     program listing number: 5.2                                        *
%*                                                                        *
%*     This program solves Newton equation for a falling object           *
%*     using Runge-Kutta 2nd and 4th order methods.                       *
%*        m = mass of the object                                          *
%*        g = acceleration due to gravity                                 *
%*        gamma = frictional coefficient                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
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
v_rk45=np.zeros(N+1)
v_ex=np.zeros(N+1)
t=

for i in range(0,N):
    
    # Runge-Kutta 2nd order
    k1 = -gamma*v_rk2[i]/m - g
    v_mid = v_rk2[i] + k1*h/2.0
    k2 = -gamma*v_mid/m - g
    v_rk2[i+1] = v_rk2[i] + k2*h
    
    # RUnge-Kutta 4th order
    k1 = -gamma*v_rk4[i]/m - g
    v_mid = v_rk4[i] + k1*h/2.0
    k2 = -gamma*v_mid/m - g
    v_mid = v_rk4[i] + k2*h/2
    k3 = -gamma*v_mid/m - g
    v_end = v_rk4[i] + k3*h    
    k4 = -gamma*v_end/m - g
    v_rk4[i+1] = v_rk4[i] + (k1+2*(k2+k3)+k4)*h/6.0
        
    # Exact solution
    v_ex[i+1] = m*g/gamma*(np.exp(-gamma*t[i+1])-1)
    
plt.ioff()
plt.figure(figsize=(12,5))

# Plot the solutions
plt.subplot(1,2,1);
plt.plot(t,v_rk2,'-b',label='RK2')
plt.plot(t,v_rk4,'-r',label='RK4')
plt.plot(t,v_ex,'-k',label='Exact')
plt.xlabel('t')
plt.ylabel('v(t)')
plt.legend(loc=1)

# Plot the absolute errors
plt.subplot(1,2,2)
plt.semilogy(t,abs(v_rk2-v_ex),'-b',label='RK2')
plt.semilogy(t,abs(v_rk4-v_ex),'-r',label='RK4')

plt.xlabel('t')
plt.ylabel('absolute error')
plt.legend(loc=3)
plt.show()