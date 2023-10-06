#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  7.1                                                       *
%*     filename: ch07pr01.py                                              *
%*     program listing number: 7.1                                        *
%*                                                                        *
%*     This program finds the first two eigen modes of standing wave in a *
%*     string using the shooting method (Numerov and secant methods).     *
%*                                                                        *
%*     Programmed by Ryoichi Kawai for Computational Physics Course.      *
%*     Last modification:  02/04/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# setting the grid
xmin=0.0
xmax=1.0
N=200
h=(xmax-xmin)/np.float(N)
x=np.linspace(xmin,xmax,N+1)
v = np.zeros(N+1)

# control parameter
kmax=100
tol=1.0e-6
delta = 1.0
found = False

lam=np.zeros(kmax+1)
err=np.zeros(kmax+1)
lam[0] = 1.0

k=0
while not(found) and k < kmax:

    # integrate ODE by Numerov method
    w = -lam[k]
    v[0]=0.0
    v[1]=delta
    for n in range(1,N):
        v[n+1] = 2.0*(1.0-5.0*h**2*w/12.0)*v[n] - (1.0+h**2*w/12.0)*v[n-1]
        v[n+1] = v[n+1]/(1.0+h**2*w/12.0)


    # error in the boundary condition
    err[k] = v[N]

    if np.abs(err[k]) < tol:
        found = True
    else:
        # secant method to guess next lambda
        if k == 0:
            lam[k+1] = lam[k]-0.1
        else:
            lam[k+1] = lam[k]-(lam[k]-lam[k-1])/(err[k]-err[k-1])*err[k]
        k+=1

s=max(v)
v=v/s
print('lambda = {0:12.7f}, exact={1:12.7f}'.format(lam[k],-np.pi**2))

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(x,v,'-b',label=r"\lambda_0=1$")
plt.xlabel('x',fontsize=14)
plt.ylabel('v(x)',fontsize=14)

plt.subplot(1,2,2)
plt.plot(np.linspace(0,k,k+1),lam[0:k+1],'-ob')
plt.plot([0,k],[-np.pi**2, -np.pi**2],'--k')
plt.plot([0,k],[-(2*np.pi)**2, -(2*np.pi)**2],'--k')
plt.xlabel('iteration',fontsize=14)
plt.ylabel(r'$\lambda$',fontsize=14)


found = False

lam[0]= -30

k=0
while not(found) and k < kmax:
    w = -lam[k]
    v[0]=0.0
    v[1]=delta
    for n in range(1,N):
        v[n+1] = 2.0*(1.0-5.0*h**2*w/12.0)*v[n] - (1.0+h**2*w/12.0)*v[n-1]
        v[n+1] = v[n+1]/(1.0+h**2*w/12.0)

    err[k] = v[N]
    if abs(err[k]) < tol:
        found = True
    else:
        if k == 0:
            lam[k+1] = lam[k]-0.1
        else:
            lam[k+1] = lam[k] -(lam[k]-lam[k-1])/(err[k]-err[k-1])*err[k]

        k+=1


s = max(v)
v = v/s
print('lambda = {0:12.7f}, exact={1:12.7f}'.format(lam[k],-(2*np.pi)**2))

plt.subplot(1,2,1)
plt.plot(x,v,'-r',label=r"\lambda_0=-30$")
plt.legend(loc=3)

plt.subplot(1,2,2)
plt.plot(np.linspace(0,k,k+1),lam[0:k+1],'-or')

plt.show()
