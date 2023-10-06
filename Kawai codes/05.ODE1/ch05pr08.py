#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 5.5.2                                                      *
%*     filename: ch05pr08.py                                              *
%*     program listing number: 5.8                                        *
%*                                                                        *
%*     This program solves the Maxwell-Bloch model of laser dynamics      *
%*     using Runge-Kutta 4th order methods.                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/22/2017.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# system parametes uncomment ONLY the desired set
gamma1_set = {0: 0.1, 1: 0.1, 2: 1.0}
gamma2_set = {0: 2.0, 1: 10.0, 2: 0.1}
gamma3_set = {0: 3.0, 1: 0.25, 2: 0.25}
c1_set = {0: 0.25, 1: 1.0, 2: 1.0}
c2_set = {0: 0.2, 1: 0.5, 2: 0.1}
c3_set = {0: 1.0, 1: 1.0, 2: 1.0}

param=np.int(input('Choose a parameter set [0-2] ='))
gamma1=gamma1_set.get(param)
gamma2=gamma2_set.get(param)
gamma3=gamma3_set.get(param)
c1=c1_set.get(param)
c2=c2_set.get(param)
c3=c3_set.get(param)

lam=np.float(input('Enter a value for lambda = '))

# Control parameters
tmax=500; N=5000; h=tmax/N
t=np.linspace(0,tmax,N+1)

E=np.zeros(N+1)
P=np.zeros(N+1)
D=np.zeros(N+1)

# initial conditions
E[0]=1.0; P[0]=1.0; D[0]=1.0

# 2nd-order Runge-Kutta method
for n in range(0,N):
    FE_n=-gamma1*E[n]+c1*P[n]
    FP_n=-gamma2*P[n]+c2*E[n]*D[n]
    FD_n=-gamma3*(D[n]-lam)-c3*E[n]*P[n]
    E_mid = E[n]+FE_n*h/2.0
    P_mid = P[n]+FP_n*h/2.0
    D_mid = D[n]+FD_n*h/2.0
    FE_mid=-gamma1*E_mid+c1*P_mid
    FP_mid=-gamma2*P_mid+c2*E_mid*D_mid
    FD_mid=-gamma3*(D_mid-lam)-c3*E_mid*P_mid
    E[n+1]=E[n]+FE_mid*h
    P[n+1]=P[n]+FP_mid*h
    D[n+1]=D[n]+FD_mid*h


# plot the dynamics of E
plt.ioff()
fig=plt.figure(figsize=(12,5))
ax=fig.add_subplot(1,2,1)
ax.plot(t,E)
ax.set_xlabel('t')
ax.set_ylabel('E(t)')

# plot 3D phase trajectory
ax = fig.add_subplot(1,2,2,projection='3d')
ax.plot(E, D, P)
ax.set_xlabel('E')
ax.set_ylabel('D')
ax.set_zlabel('P')
plt.show()
