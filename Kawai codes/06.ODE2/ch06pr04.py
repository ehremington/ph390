#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 6.3.2                                                      *
%*     filename: ch06pr04.py                                              *
%*     program listing number: 6.4-1                                      *
%*                                                                        *
%*     This program solves one-dimensional heat equation and finds        *
%*     temperature profile and heat energy transaction using              *
%*     Numerov integration.                                               *
%*     Use function: numerov_heat(TL,delta,L)                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

def numerov_heat(TL,delta,L):
    global mu
    # control parameters
    N=10000; h=L/N
    x=np.linspace(0,L,N+1)
    u=np.zeros(N+1)

    # define w(x) in Numerov method
    w=mu

    # initial conditions
    u[0]=TL

    # we guess u(h)
    u[1]=u[0]-delta

    # shoot out to x=L by the Numerov method
    for n in range(1,N):
        u[n+1] = 2.0*(1.0-5.0*h**2*w/12.0)*u[n] - (1.0+h**2*w/12.0)*u[n-1]
        u[n+1] = u[n+1]/(1+h**2*w/12.0)
        
    return [x,u]


if __name__ == "__main__":
        
    # parameters
    global mu
    mu = -10.0
    TL=10.0;  TR=0.0; L=1.0
        
    # tolerance
    tol=1.0e-9
        
    # control variable 
    found = False

    y1 = np.zeros(1000)
    y2 = np.zeros(1000)
    # first guess of delta
    y1[0] = 1
    # get the u(x)
    [x, u] = numerov_heat(TL,y1[0],L)
    N = u.size-1   # check how many grid points are used.
    y2[0]=(u[N]-TR)**2
    if abs(y2[0]) < tol:
        found = True


    if not(found):
        #second guess of delta
        y1[1] = y1[0]+0.01
        # get u(x)
        [x,u] = numerov_heat(TL,y1[1],L)
        y2[1] =(u[N]-TR)**2
        if abs(y2[1]) < tol:
            found = True

    # secant iteration
    n=1
    while not(found):
        # guess delta by secant method
        y1[n+1] = y1[n] - (y1[n]-y1[n-1])/(y2[n]-y2[n-1])*y2[n]
        # derivative of phi(x) at the end point.
        [x,u] = numerov_heat(TL,y1[n+1],L)
        y2[n+1]=(u[N]-TR)**2
        if abs(y2[n+1]) < tol:
            found = True

        n+=1


    # Energy conservation
    Q_in = -(u[1]-u[0])/(x[1]-x[0])
    Q_out= +(u[n]-u[n-1])/(x[n]-x[n-1])
    Q_diss = mu*sum(u[0:n-1:2]+4.0*u[1:n:2]+u[2:n+1:2])*(x[1]-x[0])/3.0
    print('Q_in={0:10.6f},  Q_out={1:10.6f},  Q_diss={2:10.6f},  Q_net={3:10.6f}'
          .format(Q_in, Q_out, Q_diss, Q_in+Q_out+Q_diss))

    # plot heat source
    plt.ioff()
    plt.figure(figsize=(12,5))
    plt.subplot(1,2,1)
    plt.plot(x,u,'-r')
    plt.plot([0,L],[0,0],'--k')
    plt.xlabel(r'$s$')
    plt.ylabel(r'$u(s)$')

    plt.subplot(1,2,2)
    # plot the error after each iteration.
    plt.semilogy(np.linspace(0,n,n+1),abs(y2[0:n+1]),'-o')
    plt.xlabel('iteration')
    plt.ylabel(r'$T_R$')

    plt.show()