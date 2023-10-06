#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 5.3                                                        *
%*     filename: ch05pr03.py                                              *
%*     program listing number: 5.3                                        *
%*                                                                        *
%*     This program solves Newton equation for a falling object           *
%*     using the Runge-Kutta-Fehlberg method.                             *
%*        m = mass of the object                                          *
%*        g = acceleration due to gravity                                 *
%*        gamma = frictional coefficient                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/30/2018.                                             *
%**************************************************************************
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
    

# define the right hand side of ODE   
def func(t,v):
    return -gamma*v-m*g

if __name__ == "__main__":
 
    # system parameters
    gamma=1.0
    m=1.0
    g=9.8
    # time span
    tspan=[0,10]
    # initial condition (must be ndarray)
    y0=[0]
    #relative tolerence
    rtol=1e-5
 
    # use RK45 method
    sol=solve_ivp(func,tspan,y0,method='RK45',rtol=rtol)
    # save the results
    t=sol.t
    v=list(sol.y.flat)
    
    # exact solution
    v_ex=m*g/gamma * (np.exp(-gamma*t)-1)
    
    plt.ioff()
    plt.figure(figsize=(12,5))

    plt.subplot(1,2,1) 
    plt.plot(t,v,'or',label="RK45")
    plt.plot(t,v_ex,'-k',label="Exact")
    plt.plot(t,t*0,'ob',label="Time step")
    plt.xlabel('t')
    plt.ylabel('velocity')
    plt.legend(loc=0)

    plt.subplot(1,2,2)
    plt.semilogy(t,abs(v-v_ex),'-k')
    plt.xlabel('t')
    plt.ylabel("absolute error")
    plt.show()