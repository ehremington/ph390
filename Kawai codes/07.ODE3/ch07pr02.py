#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 7.2.1                                                      *
%*     filename: ch07pr02.py                                              *
%*     program listing number: 7.2                                        *
%*                                                                        *
%*     This program finds an eigenvalue and eigenfunction of a quantum    *
%*     harmonic oscillator within a given bracket using the shooting      *
%*     method (Numerov and bisection methods).                            *
%*     Parity symmetry is taken into account.                             *
%*                                                                        *
%*     Programmed by Ryoichi Kawai for Computational Physics Course.      *
%*     Last modification:  02/04/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

def qmho_numerov(E,xmax,h):
    N=round(xmax/h)
    xmax = h*N
    psi=np.zeros(N+1)
    x =np.linspace(-xmax,0.0,N+1)
    w = -x**2+E

    psi[0]=0
    psi[1]=0.001

    # shoot out to x=L by the Numerov method
    for n in range(1,N):
        psi[n+1] = 2.0*(1.0-5.0*h**2*w[n]/12.0)*psi[n] \
            - (1.0+h**2*w[n-1]/12.0)*psi[n-1]
        psi[n+1] = psi[n+1]/(1+h**2*w[n+1]/12.0)

    return [x,psi]

if __name__ == "__main__":

    E=np.zeros(2)
    E[0]=np.float(input('Energy Lower Blacket ='))
    E[1]=np.float(input('Energy Upper Blacket ='))

    # control parameter
    symmetric = False;
    anti_symmetric = False
    found = False
    xmin=0.0
    xmax=5.0
    tol = 1e-8
    h = 0.001

    # Initial Lower bound
    [x,psi] = qmho_numerov(E[0],xmax,h)
    N=x.size-1
    error_L=(psi[N]-psi[N-1])/(x[N]-x[N-1])
    error_L2=psi[N]
    if np.abs(error_L) < tol:
        found = True
        symmetric = True
        EM = E[0]
    elif np.abs(error_L2) < tol:
        found = True
        anti_symmeric = True
        EM = E[0]

    # Initial Upper bound
    if not(found):
        [x,psi] = qmho_numerov(E[1],xmax,h)
        error_U=(psi[N]-psi[N-1])/(x[N]-x[N-1])
        error_U2=psi[N]

        if np.abs(error_U) < tol:
            found = True
            symmetric = True
            EM = E[1]
        elif np.abs(error_U2) < tol:
            found = True
            anti_symmeric = True
            EM = E[1]

        if error_U*error_L<0:
            symmetric = True
        elif error_U2*error_L2<0:
            anti_symmetric = True
            error_U = error_U2
            error_L = error_L2
        else:
            sys.exit('Blacket error!');

        if symmetric & anti_symmetric:
            sys.exit('Blacket error2!')

    # Begin bisection
    while not(found):
        EM = E.sum()*0.5
        [x,psi] = qmho_numerov(EM,xmax,h)
        if symmetric:
            error_M=(psi[N]-psi[N-1])/(x[N]-x[N-1])
        else:
            error_M=psi[N]

        if np.abs(error_M)<tol:
            found = True
        else:
            if error_M*error_L < 0:
                E[1] = EM
                error_U = error_M
            else:
                E[0] = EM
                error_L = error_M

    # output the result
    if symmetric:
        print('Symmetric state: ')
    elif anti_symmetric:
        print('Anti-Symmetric state: ')

    print('Eigenvalue= {0:12.6f}'.format(EM))

    X=np.zeros(2*N+1)
    Y=np.zeros(2*N+1)
    X[0:N+1] = x[0:N+1]
    X[N+1:2*N]=-x[N-1:0:-1]
    if symmetric:
        Y[0:N+1] = psi[0:N+1];
        Y[N+1:2*N]=psi[N-1:0:-1]
    else:
        Y[1:N] = psi[1:N]
        Y[N+1:2*N]=-psi[N-1:0:-1]

    A = sum(Y[0:2*N-3:2]**2+4.0*Y[1:2*N-2:2]**2+Y[2:2*N-1:2]**2)*h/3.0
    Y = Y / np.sqrt(A)
    plt.figure(figsize=(6,5))
    plt.plot(X,Y,'-b')
    plt.plot([-xmax,xmax],[0,0],'-k')
    plt.plot([0,0],[-1,1],'-k')
    plt.xlabel('x')
    plt.ylabel(r'$\psi(x)$')
    plt.show()
