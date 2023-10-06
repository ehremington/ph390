#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  6.1                                                       *
%*     filename: ch06pr01.py                                              *
%*     program listing number: 6.1-1                                      *
%*                                                                        *
%*     This program determines a launching speed that a rocket necessary  *
%*     to reach height yf in travel time tf.                              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

def f(v):
    # right hand side of the ODE
    g=9.8; m=1.0; C=0.01
    return -(C/m)*np.abs(v)*v-g

def rocket_trajectory(vi,t):
    # Solve the ODE using RK5 and return the final position abd velocity
    N=1000
    h=t/N
    y0=0.0
    v0=vi
    
    for n in range(N):
        ky1 = v0
        kv1 = f(v0)

        y_mid = y0 + ky1*h/2.0
        v_mid = v0 + kv1*h/2.0
        ky2 = v_mid
        kv2 = f(v_mid)
      
        y_mid = y0 + ky2*h/2.0
        v_mid = v0 + kv2*h/2.0
        ky3 = v_mid
        kv3 = f(v_mid)

        y_end = y0 + ky3*h
        v_end = v0 + kv3*h
        ky4 = v_end
        kv4 = f(v_end)
    
        y0=y0+(ky1+2*(ky2+ky3)+ky4)*h/6.0
        v0=v0+(kv1+2*(kv2+kv3)+kv4)*h/6.0

    return [y0,v0]

if __name__ == "__main__":
    # set the boundary conditions
    yf=100.0; tf=2.0

    # tolerance
    tol=1e-8

    # control variable 
    nmax = 100
    found = False

    y=np.zeros(nmax+1)
    v=np.zeros(nmax+1)
    # first guess
    n=1
    v[n] = 50.0
    [y[n], vf] = rocket_trajectory(v[n],tf)
    if np.abs(y[n]-yf) < tol :
        found = True
        v0 = v[n]

    #second guess
    n+=1
    v[n] = 51.0
    [y[n], vf] = rocket_trajectory(v[n],tf)
    if np.abs(y[n]-yf) < tol:
        found = True
        v0 = v[n]

# secant iteration
    while not(found) :
        v[n+1] = v[n] - (v[n]-v[n-1])/(y[n]-y[n-1])*(y[n]-yf)
        [y[n+1], vf] = rocket_trajectory(v[n+1],tf)
        if np.abs(y[n+1]-yf) < tol:
            found = True
            v0 = v[n+1]
        n+=1

    # show the result
    print('initial velocity = {0:10.6f} final velocity = {1:10.6f}'
          .format(v[n],vf))

    # plot the convergency
    plt.ioff()
    plt.figure(figsize=(6,5))
    plt.plot(np.linspace(1,n,n),v[1:n+1],'-ob',label='numerical')
    plt.plot([0,n+2],[101.9281,101.9281],'--',label='exact')
    plt.xlabel('Iteration')
    plt.ylabel('$v_0$')
    plt.legend(loc=4)
    plt.show()