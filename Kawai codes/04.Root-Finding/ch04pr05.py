#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*  Example 4.5                                                           *
%*  filename: ch04pr05.py                                                 *
%*  program listing number: 4.5                                           *
%*                                                                        *
%*     This program calculates magnetization as a function of             *
%*     temperature using the mean field Ising model:                      *
%*         x = x tan(x/S)                                                 *
%*     where the variables are normalized as                              *
%*         x = m/m_0 and S=k*T/(C*m_0) .                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/14/2017.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

def f(x,a):
    return x - np.tanh(a*x)

def df(x,a):
    return 1.0 - a/np.cosh(a*x)**2

if __name__ == "__main__": 
    N=100 # number of data points
    tol=1.0e-6  # tolerance

    m=np.zeros(N+1)
    t=np.zeros(N+1)
    dS = 2.0/N   # step size in temperature

    
    for k in range(0,N):
        S=(k+1)*dS  # temperature
        a=1.0/S
        # Newton-Raphson method
        x=2.0    # initial guess
        fx=f(x,a)
        while fx>tol:
            x=x-fx/df(x,a)
            fx=f(x,a)

        # store the magnetization and temperature  
        m[k]=x
        t[k]=S

    plt.ioff()
    plt.figure(figsize=(12,5))
    plt.plot(t,m, '-r')
    plt.plot(t,-m, '-r')
    plt.plot([0,1],[0,0],'--w')
    plt.ylim([-1.2,1.2])
    plt.xlabel('T')
    plt.ylabel('m')
    plt.savefig('ch04r05.pdf')
    plt.show()
 