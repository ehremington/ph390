#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 16:19:49 2014

@author: eric
"""

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
        x[n+1:] = [x[n] + h/6*(k1+2*(k2+k3)+k4)]
        t[n+1:] = [t[n] + h]
        n = n + 1
    return x,t


def coupledRK2(f,g,h,ending_time,x0 = 0,y0 = 0, t0 = 0):
    x = [x0]
    y = [y0]
    t = [t0]
    n = 0
    while t[n]<ending_time:
        k1 = f(x[n],y[n],t[n])
        j1 = g(x[n],y[n],t[n])
        k2 = f(x[n]+k1*h/2,y[n]+j1*h/2,t[n]+h/2)
        j2 = g(x[n]+k1*h/2,y[n]+j1*h/2,t[n]+h/2)
        x[n+1:] = [x[n] + k2*h]
        y[n+1:] = [y[n] + j2*h]
        t[n+1:] = [t[n] + h]
        n = n + 1
    return x,y,t    

def coupledRK4(f,g,h,ending_time,x0 = 0,y0 = 0, t0 = 0):
    x = [x0]
    y = [y0]
    t = [t0]
    n = 0
    while t[n]<ending_time:
        k1 = f(x[n],y[n],t[n])
        j1 = g(x[n],y[n],t[n])
        k2 = f(x[n]+k1*h/2,y[n]+j1*h/2,t[n]+h/2)
        j2 = g(x[n]+k1*h/2,y[n]+j1*h/2,t[n]+h/2)
        k3 = f(x[n]+k2*h/2,y[n]+j2*h/2,t[n]+h/2)
        j3 = g(x[n]+k2*h/2,y[n]+j2*h/2,t[n]+h/2)
        k4 = f(x[n]+k3*h,y[n]+j3*h,t[n]+h)
        j4 = g(x[n]+k3*h,y[n]+j3*h,t[n]+h)
        x[n+1:] = [x[n] + h/6*(k1+2*(k2+k3)+k4)]
        y[n+1:] = [y[n] + h/6*(j1+2*(j2+j3)+j4)]
        t[n+1:] = [t[n]+h]
        n = n + 1
    return x,y,t
    
    
    
    
    
    
    
    
    