#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 7.2.3                                                      *
%*     filename: ch07pr04.py                                              *
%*     program listing number: 7.4                                        *
%*                                                                        *
%*     This program finds an eigenvalue and eigenfunction of a quantum    *
%*     particle in the Morse potential using the shooting.                *
%*     method (Numerov and secant methods).                               *
%*                                                                        *
%*     Programmed by Ryoichi Kawai for Computational Physics Course.      *
%*     Last modification:  02/04/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt
import sys


# Use the atomic unit
mp=1836.4
e=27.2114
a0=0.529177
# parameters for H2
D=4.75/e; R0=0.742/a0; a=1.44/R0; name='H_2'; m=mp; mu=m/2.

# parameters for I2
#D=1.56/e; R0=2.66/a0; a=4.94/R0; name='I_2'; m=127*mp; mu=m/2;

# zero point energy
omega=np.sqrt(2.0*D*a**2/mu)
DE=omega/5.0
d=1.0/np.sqrt(mu*omega)*12.0


# define the discrete coordinate
N=1000
h=d/N
x=np.linspace(-d/2.0,d/2.0,N+1)
psi=np.zeros(N+1)
phi=np.zeros((3,N+1))
# evaluate the potential
U=D*(np.exp(-2.0*a*x)-2.0*np.exp(-a*x))


eigval=np.zeros(3)
E_M=-D

for L in [0, 1, 2]:

    # Initial Bracketting
    E_L = E_M+DE
    w = 2*mu*(E_L-U)
    psi[0]=0.0
    psi[1]=0.1
    # shoot out to x=L by the Numerov method
    for n in range(1,N):
        psi[n+1] = 2.0*(1.0-5.0*h**2*w[n]/12.0)*psi[n] \
            - (1.0+h**2*w[n-1]/12.0)*psi[n-1]
        psi[n+1] = psi[n+1]/(1+h**2*w[n+1]/12.0)

    ERR_L=psi[N]
    found=False

    while not(found):

        E_U = E_L+DE
        if E_U>0.0:
            sys.exit('Blacket error')

        w = 2.0*mu*(E_U-U)
        psi[0]=0.0
        psi[1]=0.1
        for n in range(1,N):
            psi[n+1] = 2.0*(1.0-5.0*h**2*w[n]/12.0)*psi[n] \
                - (1.0+h**2*w[n-1]/12.0)*psi[n-1]
            psi[n+1] = psi[n+1]/(1+h**2*w[n+1]/12.0)

        ERR_U=psi[N]
        if np.sign(ERR_L)*np.sign(ERR_U)<0.0:
            found=True
        else:
            E_L=E_U
            ERR_L=ERR_U

    found=False
    tol1=1e-8
    tol2=1e-15

    while not(found):
        E_M=(E_L+E_U)/2.0
        w = 2.0*mu*(E_M-U)
        psi[0]=0.0
        psi[1]=0.1

        for n in range(1,N):
            psi[n+1] = 2.0*(1.0-5.0*h**2*w[n]/12.0)*psi[n] \
                - (1.0+h**2*w[n-1]/12.0)*psi[n-1]
            psi[n+1] = psi[n+1]/(1+h**2*w[n+1]/12.0)

        ERR_M=psi[N]
        if np.abs(ERR_M) < tol1 or np.abs(E_L-E_M)<tol2:
            found = True

        if np.sign(ERR_M)*np.sign(ERR_L) > 0.0 :
            E_L=E_M
        else:
            E_U=E_M

    eigval[L] = E_M
    # normalize the wave function
    z=sum(psi**2)*h
    phi[L,:]=psi/np.sqrt(z);
    exact=-D*(1.0-a/np.sqrt(2.0*mu*D)*(L+1./2.))**2
    print('n={0:3d}, E={1:10.7f}, Exact={2:10.7f}'.format(L-1,E_M,exact))


# Ploting potential and eigenvalues
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
N0=np.int(np.ceil(N/5))
plt.plot([x[N0],x[N]],[eigval[0],eigval[0]],'-b',label='n=0')
plt.plot([x[N0],x[N]],[eigval[1],eigval[1]],'-g',label='n=1')
plt.plot([x[N0],x[N]],[eigval[2],eigval[2]],'-r',label='n=2')
plt.plot(x[N0:N],U[N0:N],'-k')
plt.xlabel(r'$\xi$',fontsize=14)
plt.ylabel(r'$U(\xi)$',fontsize=14)
plt.legend(loc=1)

# Ploting wave functions
plt.subplot(1,2,2)
plt.plot(x,phi[0,:],'-b')
plt.plot(x,phi[1,:],'-g')
plt.plot(x,phi[2,:],'-r')
plt.xlabel(r'$\xi$',fontsize=14)
plt.ylabel(r'$\psi(\xi)$',fontsize=14)

plt.show()
