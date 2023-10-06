#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 18.4                                                       *
%*     filename: ch18pr04.py                                              *
%*     program listing number: 18.4                                       *
%*                                                                        *
%*     This program simulates the Orstein-Uhlenbeck process.              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/26/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# system parameters
gamma=0.1    # friction
T=1.0        # temperature
k=1.0        # spring constant

# control parameter
tau=1000.0
dt=0.01
ds=sqrt(2*T/gamma*dt)
dk=k*dt
N=np.int(tau/dt)   # numner of time steps
M=1000             # number of realization

x=np.zeros(N+1)
t=np.linspace(0.0,dt*N,N+1)

sumx1=zeros(N+1)
sumx2=zeros(N+1)

for j in range(0,M):
    
    x[0]=0                # initial position
    g=np.random.randn(N)  # Gaussian random numner

    for i in range(0,N):
        x[i+1] = x[i]+ds*g[i] - dk*x[i]               # Eulrer step
        x[i+1] = x[i]+ds*g[i] - dk*(x[i+1]+x[i])/2.0  # Huen step

    sumx1 = sumx1+x
    sumx2 = sumx2+x**2

plt.close('all')
plt.figure(figsize=(6,5))
plt.plot(t,x,'-b')
plt.plot([t[0],t[N]],[0,0],'--r')
plt.xlabel('t',fontsize=14)
plt.ylabel('x(t)',fontsize=14)
plt.show()

plt.figure(figsize=(6,5))
xmean=sumx1/M
sumx2=sumx2/M
sigma2=sumx2-xmean**2
theory = T/gamma
plt.plot(t,xmean,'-k',label=r'$\langle x \rangle$')
plt.plot(t,sigma2,'-b',label=r'$\sigma^2$')
plt.plot([t[0],t[N]],[theory,theory],'--r',label='exact')
plt.xlabel('t',fontsize=14)
plt.ylabel('moments',fontsize=14)
plt.legend(loc=5)
plt.show()

plt.figure(figsize=(6,5))
counts, bins = np.histogram(x,21,normed=1)
X=(bins[:-1] + bins[1:])/2.0
Y=sqrt(gamma/(2.0*np.pi*T))*exp(-k*X**2*gamma/(2*T))
plt.plot(X,counts,'ok',label='simulation')
plt.plot(X,Y,'-r',label='exact')
plt.xlabel('x',fontsize=14)
plt.ylabel(r'$\rho(x)$',fontsize=14)
plt.show()
