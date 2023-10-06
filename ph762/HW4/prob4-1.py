#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 16:19:49 2014

@author: eric
"""
import matplotlib.pyplot as plt
#import math
import numpy as np

#define initial conditions
g = 0.1
m = 1
v_0 = 10
tau = 2

t_0 = 0

#the following has explicit x and t placements
def f(v,t):
    return -g/m*v*np.exp(-t/tau)

#euler method ode solver
def euler(function, h, ending_time, t0 = 0, x0 = 1):
    x = [x0]  #initializing lists and setting initial conditions
    t = [t0]    
    n=0
    while t[n]<ending_time:
        x[n+1:] = [x[n] + function(x[n],t[n])*h]
        t[n+1:] = [t[n]+h]
        n = n+1
    return x,t

def RungeKutta2nd(function,h,ending_time,t0 = 0, x0 = 1):
    x = [x0]
    t = [t0]
    n = 0
    while t[n]<ending_time:
        k1 = function(x[n],t[n])
        k2 = function(x[n]+k1*h/2,t[n]+h/2)
        x[n+1:] = [x[n]+k2*h]
        t[n+1:] = [t[n]+h]
        n = n+1
        
    return x,t
    
    
def RungeKutta4th(function,h,ending_time,t0 = 0, x0 = 1):
    x = [x0]
    t = [t0]
    n = 0
    while t[n]<ending_time:
        k1 = function(x[n],t[n])
        k2 = function(x[n]+k1*h/2,t[n]+h/2)
        k3 = function(x[n]+k2*h/2,t[n]+h/2)
        k4 = function(x[n]+k3*h,t[n]+h)
        x[n+1:] = [x[n] + h/6*(k1+2*k2+2*k3+k4)]
        t[n+1:] = [t[n] + h]
        n = n + 1
    return x,t
        
x4th,t4th = np.asarray(RungeKutta4th(f,.1,20,0,v_0))

x2nd,t2nd = np.asarray(RungeKutta2nd(f,.1,20,0,v_0))  #performs the euler method and converts output to numpy array

xe = v_0*np.exp(-g*tau/m*(1-np.exp(-t4th/tau))) #exact solution determined analytically

diff4th = abs(xe - x4th)
diff2nd = abs(xe - x2nd)

plt.close()  #clears any current figure
f, ax = plt.subplots(2, sharex = True)  #shared x axis
ax[0].plot(t4th,x4th,'o',t2nd,x2nd,'v')
ax[0].set_ylabel('velocity (m/s)')
ax[0].set_title('2nd and 4th order Runge-Kutta')
ax[1].plot(t4th,diff4th,'o',t2nd,diff2nd,'v')
ax[1].set_title('Error Comparison')
ax[1].set_ylabel('Error from exact solution (m/s)')
plt.xlabel('time (s)')
plt.show()

