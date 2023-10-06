#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 7.2.2                                                      *
%*     filename: ch07pr03.py                                              *
%*     program listing number: 7.3                                        *
%*                                                                        *
%*     This program finds an eigenvalue and eigenfunction of a quantum    *
%*     bouncing ball within a given bracket using the shooting            *
%*     method (Numerov and secant methods).                               *
%*                                                                        *
%*     Programmed by Ryoichi Kawai for Computational Physics Course.      *
%*     Last modification:  02/04/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# setting the grid
L=10.0; N=500; h=L/N
x=np.linspace(0,L,N+1)
psi=np.zeros(N+1)

# control parameter
tol=1e-6
delta = 0.01
found = False
kmax=100
eigval=np.zeros(kmax+1)
err=np.zeros(kmax+1)
eigval[0] = np.float(input('Initial Guess ='))  # initial guess

k=0
# secant iteration
while not(found) and k<kmax:
    # Numerov integration
    w = eigval[k]-x
    psi[N]=0.0
    psi[N-1]=delta
    for n in range(N-1,0,-1):
        psi[n-1]=2.0*(1.0-5.0*h**2*w[n]/12.0)*psi[n] \
            - (1.0+h**2*w[n+1]/12.)*psi[n+1]
        psi[n-1] = psi[n-1]/(1.0+h**2*w[n-1]/12.0)

    err[k] = psi[0]
    if np.abs(err[k]) < tol:
        found = True
    else:
        if k == 0:
            eigval[k+1] = eigval[k]-0.1  # second guess
        else:
            # suggestion by the secant method
            eigval[k+1] = eigval[k] \
                -(eigval[k]-eigval[k-1])/(err[k]-err[k-1])*err[k]

        k+=1

# normalize the solution at the maximum
psi=psi/max(psi)
print('Eigenvalue = {0:12.7f}'.format(eigval[k]))

# plot eigenfunction
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(x,psi,'-b',linewidth=2)
plt.plot([x[0],x[N]],[0,0],'-k')
plt.xlabel('x',fontsize=14)
plt.ylabel('v(x)',fontsize=14)

# plot eigenvalue
plt.subplot(1,2,2)
plt.plot(np.linspace(0,k,k+1),eigval[0:k+1],'-ob')
plt.plot([0,k],[-np.pi**2, -np.pi**2],'--k');
plt.xlabel('iteration',fontsize=14)
plt.ylabel('E',fontsize=14)

plt.show()
