#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Exercise 18.3                                                      *
%*     filename: ch18pr03.py                                              *
%*     program listing number: 18.3                                       *
%*                                                                        *
%*     This program calculates temporal autocorrelation of velocity of    *
%*     one-dimensional Brownian particles.                                *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/26/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# system parameters
m=1.0        # mass
gamma=0.1/m  # friction
T=2.0        # temperature
w=np.sqrt(2.0*gamma*T)/m

# control parameters
N=2**20  #   number of data points
dt=0.01
ds=np.sqrt(dt)*w

x=np.zeros(N)
v=np.zeros(N)
t=np.linspace(0.0,dt*(N-1),N)

# initial conditions
x[0]=0.0   # initial position
v[0]=1.0   # initial velocity

# random number
g=np.random.randn(N)

# solve Langevin Eq.
for i in range(0,N-1):
    # Euler step
    x[i+1]=x[i]+v[i]*dt;
    v[i+1]=v[i]-gamma*v[i]*dt+g[i]*ds
    #Heun step
    x[i+1]=x[i]+(v[i]+v[i+1])*dt/2.0
    v[i+1]=v[i]-gamma*(v[i]+v[i+1])*dt/2.0+g[i]*ds

# velocity autocorrelation
z=np.fft.fft(v)
y=abs(z)**2
vc=np.fft.ifft(y)/N
vc=vc/vc[0]
vx=np.exp(-gamma*t)

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(t,v,'-b')
plt.plot([0,t[N-1]],[0,0],'--r')
plt.xlabel('t',fontsize=14)
plt.ylabel('v(t)',fontsize=14)
plt.show()

plt.subplot(1,2,2)
nc=np.int(50/dt)
plt.plot(t[0:nc],vx[0:nc],'-r',label='simulation')
plt.plot(t[0:nc],vc[0:nc],'-b',label='exact')
plt.legend(loc=1)
plt.xlabel(r'$\tau$',fontsize=14)
plt.ylabel(r'$\langle v(\tau) v(0) \rangle$',fontsize=14)
plt.plot([0,t[nc]],[0,0],'--k')
plt.show()
