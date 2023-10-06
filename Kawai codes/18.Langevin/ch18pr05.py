#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 18.2.1                                                     *
%*     filename: ch18pr05.py                                              *
%*     program listing number: 18.5                                       *
%*                                                                        *
%*     This program simulates a flashing ratchet.                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/26/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# system parameters
gamma=0.1     # friction
T=1.0         # temperature
L=1.0         # periodicity
k=2.0*np.pi/L # wave number
U0=10.0       # potential strength
F0=U0*np.pi/L # force strength
Fx=2.0          # external load
h=10.0        # pontential on-off period
D=T/gamma     # noise strength

# control parameter
tau=100.       # total time
dt=0.005      # time step
ds=np.sqrt(2.0*D*dt)   # Wiener time step
N=np.int(tau/dt)    # number of iterations
M=1000              # sample numbers

# prepare the arrays
sumx1=np.zeros(N+1)
sumx2=np.zeros(N+1)
x=np.zeros(N+1)
t=np.linspace(0.0,dt*N,N+1)

plt.close('all')
plt.figure(figsize=(6,5))
for j in range(0,M):
    u=0.0
    potential=True 
    x[0]=0.0 # initial position
    g=np.random.randn(N)  # generate Gaussian white noise
    
    for i in range(0,N):
        u=u+dt
        if potential:  # diffusion inside the potential
            f1 = -F0*(2.0*np.cos(k*x[i])-np.cos(2.0*k*x[i]))+Fx
            x[i+1] = x[i]+g[i]*ds + f1*dt
            f2 = -F0*(2.0*np.cos(k*x[i+1])-np.cos(2.0*k*x[i+1]))+Fx     
            x[i+1] = x[i]+g[i]*ds + (f1+f2)*dt/2.0
        else:  # free diffusion
            x[i+1] = x[i]+g[i]*ds+Fx*dt
   
        if u>h:
            potential = not(potential)
            u=0.0

    if j < 11:  # show only 10 trajectories
        plt.plot(t,x,'-k')
        plt.pause(0.00001)

    sumx1 = sumx1+x
    sumx2 = sumx2+x**2
    
# statstical analysis
xmean=sumx1/M
sumx2=sumx2/M
sigma2=sumx2-xmean**2

plt.plot([t[0],t[N]],[0,0],'--r')
plt.plot(t,xmean,'-r')
plt.xlabel('t',fontsize=14)
plt.ylabel('x(t)',fontsize=14)
plt.show()

plt.figure(figsize=(6,5))
plt.plot(t,xmean,'-b',label=r'$\langle x \rangle$')
plt.plot(t,sigma2,'-r',label=r'$\sigma^2$')
plt.xlabel('t',fontsize=14)
plt.ylabel('moments',fontsize=14)
plt.legend(loc=5)
plt.plot([t[0],t[N]],[0,0],'--k');
plt.show()