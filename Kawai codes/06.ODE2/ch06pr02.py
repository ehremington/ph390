#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  6.2                                                       *
%*     filename: ch06pr02.m                                               *
%*     program listing number: 6.2-1                                      *
%*                                                                        *
%*     This program solves one-dimensional Poisson equation using         *
%*     Numerov integration and secant root finding methods.               *
%*     Use function: numerov_poisson(y,L)                                 *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

def S(x):
    return -2.0*x*np.exp(-x**2)
    
def numerov_poisson(y1,L,N):
#    control parameters
    h=L/N
    x=np.linspace(0,L,N+1)
    y=np.zeros(N+1)  # field phi(x)
    s=np.zeros(N+1)

    # initial conditions
    # due to symmetry phi(0)=0
    y[0]=0.0
    s[0]=S(y[0])

    # we guess phi(h)=phi_1
    y[1]=y1
    s[1]=S(x[1])

    # shoot out to x=L by the Numerov method
    n=1
    while n < N :
        s[n+1]=S(x[n+1])
        y[n+1] = 2.0*y[n]-y[n-1]+(s[n+1]+10.0*s[n]+s[n-1])*h**2/12.0
        n+=1
        
    return x, y

if __name__ == "__main__":
    # set the boundary conditions
    L=10.0
    N=10000
    # tolerance
    tol=1.0e-16
    # control variable 
    found = False

    y1=np.zeros(101)
    y2=np.zeros(101)
    n=1
    # first guess of phi_1
    y1[0] = 0.1

    # get the potential phi(x)
    x, y = numerov_poisson(y1[0],L,N)

    # derivative of phi(x) at the end point.
    y2[0] = (y[N]-y[N-1])/(x[N]-x[N-1])
    if np.abs(y2[0]) < tol:
        found = True

    if not(found):
        # second guess of phi_1
        y1[1] = y1[0]+0.01
        # get the potential phi(x)
        x, y = numerov_poisson(y1[1],L,N)
        # derivative of phi(x) at the end point.
        y2[1] = (y[N]-y[N-1])/(x[N]-x[N-1])
        if np.abs(y2[1]) < tol:
                  found = True

    # secant iteration
    n=1
    while not(found):

        # guess phi_1 by secant method
        y1[n+1] = y1[n] - (y1[n]-y1[n-1])/(y2[n]-y2[n-1])*y2[n]
        # derivative of phi(x) at the end point.
        x, y = numerov_poisson(y1[n+1],L,N)
        # derivative of phi(x) at the end point.

        y2[n+1] = (y[N]-y[N-1])/(x[N]-x[N-1])
        if np.abs(y2[n+1]) < tol:
            found = True

        n+=1
        print("Itertation ={0:5d}, y2={1:15.5e}".format(n,y2[n]))

    # construct the whole curve from x=-L to x=L.
    X=np.zeros(2*N+1)
    Y=np.zeros(2*N+1)
    X[0:N] = -x[N:0:-1]; X[N:2*N+1]=x[0:N+1]
    Y[0:N] = -y[N:0:-1]; Y[N:2*N+1]=y[0:N+1]

    #plot charge density
    plt.ioff()
    plt.figure(figsize=(12,5))
    plt.subplot(1,2,1)
    plt.plot(X,2*X*np.exp(-X*X),'-g',label=r"$\rho(x)$")

    # plot the numerical potential
    plt.plot(X,Y,'-r',linewidth=2.0,label="Numerical")
    # plot the analytic potential
    plt.plot(X,np.sqrt(np.pi)/2.0*erf(X),'-b',label="Exact")
    plt.plot([-L,L],[0.0,0.0],'--k')
    plt.xlabel('x')
    plt.ylabel(r"$\phi(x)$")
    plt.legend(loc=4)

    plt.subplot(1,2,2)
    # plot the improvment of the first point.
    plt.semilogy(np.linspace(0,n,n+1),abs(y2[0:n+1]),'-o')
    plt.xlabel('iteration')
    plt.ylabel("$\phi'(L)$")
    plt.show()
