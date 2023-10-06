#!/usr/bin/python3
# -*- coding: utf-8 -*-

import rootFinding as rf


#euler method ode solver
def euler(function, h, tf, t0 = 0, x0 = 1):
    x = [x0]  #initializing lists and setting initial conditions
    t = [t0]    
    n=0
    while t[n]<tf:
        x[n+1:] = [x[n] + function(x[n],t[n])*h]
        t[n+1:] = [t[n]+h]
        n = n+1
    return x,t

def RK2(function,h,tf,t0 = 0, x0 = 1):
    x = [x0]
    t = [t0]
    n = 0
    while t[n]<tf:
        k1 = function(x[n],t[n])
        k2 = function(x[n]+k1*h/2,t[n]+h/2)
        x[n+1:] = [x[n]+k2*h]
        t[n+1:] = [t[n]+h]
        n = n+1
        
    return x,t
    
    
def RK4(function,h,tf,t0 = 0, x0 = 1):
    x = [x0]
    t = [t0]
    n = 0
    while t[n]<tf:
        k1 = function(x[n],t[n])
        k2 = function(x[n]+k1*h/2,t[n]+h/2)
        k3 = function(x[n]+k2*h/2,t[n]+h/2)
        k4 = function(x[n]+k3*h,t[n]+h)
        x[n+1:] = [x[n] + h/6*(k1+2*(k2+k3)+k4)]
        t[n+1:] = [t[n] + h]
        n = n + 1
    return x,t


def cRK2(f,g,tf,y0 = 0,x0 = 0, t0 = 0,h=.01):
    x = [x0]
    y = [y0]
    t = [t0]
    n = 0
    while t[n]<tf:
        k1 = f(x[n],y[n],t[n])
        j1 = g(x[n],y[n],t[n])
        k2 = f(x[n]+k1*h/2,y[n]+j1*h/2,t[n]+h/2)
        j2 = g(x[n]+k1*h/2,y[n]+j1*h/2,t[n]+h/2)
        x[n+1:] = [x[n] + k2*h]
        y[n+1:] = [y[n] + j2*h]
        t[n+1:] = [t[n] + h]
        n = n + 1
    return x,y,t    

def cRK4(f,g,tf,y0 = 0,x0 = 0, t0 = 0,h=2**-10):
    x = [x0]
    y = [y0]
    t = [t0]
    n = 0
    while t[n]<tf:
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
        n = n + 1
        t[n+1:] = [t[0]+n*h]     
    return x,y,t


def cRK4core(f,g,y0=0,x0=0,t0=0,h=2**-10):
    x = x0
    y = y0
    t = t0
    k1 = f(x,y,t)
    j1 = g(x,y,t)
    k2 = f(x+k1*h/2,y+j1*h/2,t+h/2)
    j2 = g(x+k1*h/2,y+j1*h/2,t+h/2)
    k3 = f(x+k2*h/2,y+j2*h/2,t+h/2)
    j3 = g(x+k2*h/2,y+j2*h/2,t+h/2)
    k4 = f(x+k3*h,y+j3*h,t+h)
    j4 = g(x+k3*h,y+j3*h,t+h)
    x = x + h/6*(k1+2*(k2+k3)+k4)
    y = y + h/6*(j1+2*(j2+j3)+j4)
    t = t+h
    return x,y,t
    
def cRK4finalpoint(f,g,tf,y0=0,x0=0,t0=0,h=2**-10):
    x = x0
    y = y0
    t = t0
    while t<tf:
        k1 = f(x,y,t)
        j1 = g(x,y,t)
        k2 = f(x+k1*h/2,y+j1*h/2,t+h/2)
        j2 = g(x+k1*h/2,y+j1*h/2,t+h/2)
        k3 = f(x+k2*h/2,y+j2*h/2,t+h/2)
        j3 = g(x+k2*h/2,y+j2*h/2,t+h/2)
        k4 = f(x+k3*h,y+j3*h,t+h)
        j4 = g(x+k3*h,y+j3*h,t+h)
        x = x + h/6*(k1+2*(k2+k3)+k4)
        y = y + h/6*(j1+2*(j2+j3)+j4)
        t = t+h
    return x,y,t
    
def numerov(s,w,tf,y0=0,x0=0,t0=0,h=2**-12):
    def f(x,y,t):
        return y
    def g(x,y,t):
        return s(t)-w(t)*x
    
    """
    k1 = y0
    k2 = y0 + h/2*(s(t0)-w(t0)*x0)
    k3 = (1+h/2)*y0 + h/2*(s(t0+h/2)-w(t0+h/2)*(x0+y0*h/2))
    k4 = y0 + h*(s(t0+h/2)-w(t0+h/2)*((1-h/2)*x0-h*h/4*(s(t0)-w(t0)*x0)))"""
    t = [t0,t0+h]
    w_list = [w(t[0]),w(t[1])]
    s_list = [s(t[0]),s(t[1])]
    x1,y1,t1 = cRK4core(f,g,y0,x0,t0,h)
    #x = [x0,x1]
    #x = [x0,x0+h/6*(k1+2*(k2+k3)+k4)]
    #x = [x0,x0 + h/6*(y0*(4+2*h+h*h/2)+(s_list[0]-w_list[0]*x0)*(h+h*h/2+h*h*h/4))]
    x = [x0,x0+y0*h]
    #x = [x0,x0+y0*h+.5*h*h*(s(t0)-w(t0)*x0)]
    #x = [x0,np.cos(h)]
    n = 1
    
    
    while t[n]<tf:
        t[n+1:] = [t[1]+n*h]
        w_list[n+1:] = [w(t[n+1])]
        s_list[n+1:] = [s(t[n+1])]
        x[n+1:] = [1/(1+h*h/12*w_list[n+1])*\
                   (2*(1-5/12*h*h*w_list[n])*x[n]-(1+h*h/12*w_list[n-1])*x[n-1]+h*h/12*(s_list[n+1]+10*s_list[n]+s_list[n-1]))]
        
        n = n+1
    return x,t

    
def shooting(F, t_f, x_f, t_0, x_0, y0_guess1, y0_guess2, tolerance = 2**-32):
    """
    the initial condition list must be in the order of tf,xf,t0,x0
    """   
    def f(x,y,t): return y
    def g(x,y,t): return F(x,y,t)
    def shot(y):
        xf_test, yf_test, tf_test = cRK4finalpoint(f,g,t_f,y,x_0,t_0)
        accuracy = x_f-xf_test
        return accuracy
    y_0_temp,n_temp,bracket = rf.bisection(shot,y0_guess1,y0_guess2,2**-3)
    y_0,n = rf.secant(shot,y_0_temp,(bracket[1]-bracket[0])/10,tolerance)
    
    return y_0