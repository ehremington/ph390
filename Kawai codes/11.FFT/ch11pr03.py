#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 11.4.3                                                     *
%*     filename: ch11pr03.m                                               *
%*     program listing number: 11.3                                       *
%*                                                                        *
%*     This program calculates the power spectrum of coupled oscillators. *
%*                                                                        *
%*     Uses MATLAB function fft()                                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2015.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# system parameters
m=[2.,4.,3]
k=[2.,4.,4.,2.]
L=[1.,2.,1.,2.,8.];
K=[[k[0]+k[1],-k[1],0], [-k[1],k[1]+k[2],-k[2]], [0,-k[2],k[2]+k[3]]]
b=[k[0]*L[0]-k[1]*L[1], k[1]*L[1]-k[2]*L[2],k[2]*L[2]+k[3]*(L[4]-L[3])]
M=[[1./m[0],0.,0.],[0.,1./m[1],0.],[0.,0.,1./m[2]]]
b=np.array(b)
K=np.array(K)
M=np.array(M)

# initial conditions
x0=np.array([5./3.,4.,16./3.])   # at the equilibrium
v0=np.array([1.,0.,0.])

# control parameters
T=200.
N=4096
dt=T/N
dw=2.*np.pi/T
w=np.linspace(0.,dw*(N-1),N)
t=np.linspace(0.,dt*(N-1),N)
x=np.zeros((3,N))
v=np.zeros((3,N))
f=np.zeros(3)
X=np.zeros(N)
Y=np.zeros(N)
# First we calculate the trajectory of oscillators
x[:,0]=x0
v[:,0]=v0

#1st step (Euler step)
x[:,1] = x[:,0]+v[:,0]*dt

# Verlet algorithm
for i in range(2,N):
    f= -np.dot(K,x[:,i-1])+b
    x[:,i]=2.0*x[:,i-1]-x[:,i-2]+np.dot(M,f)*dt**2
    v[:,i-1]=(x[:,i]-x[:,i-2])/(2.0*dt)

# Now, we analyze the trajectories.
X=x[0,:]-x0[0] # deviation form the equilibrium
Y=np.fft.fft(X)/T  # Fourer transform
S = np.abs(Y**2) # power spectrum

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(t,x[0,:],'-b',label=r'$x_1$')
plt.plot(t,x[1,:],'-g',label=r'$x_2$')
plt.plot(t,x[2,:],'-r',label=r'$x_3$')
plt.xlabel(r'$t$',fontsize=14)
plt.ylabel(r'$x(t)$',fontsize=14)
plt.axis([0, T, 0, 7,])
plt.legend(loc=2)

plt.subplot(1,2,2)
plt.plot(w,S,'-b')
plt.axis([0, 3, 0, 8,])
plt.xlabel(r'$\omega$',fontsize=14)
plt.ylabel(r'$S(\omega)$',fontsize=14)
plt.plot( [2.056127, 2.056127], [0, 8], '--k')
plt.plot( [1.540700, 1.540700], [0, 8], '--k')
plt.plot( [0.631338, 0.631338], [0, 8], '--k')
plt.show()
