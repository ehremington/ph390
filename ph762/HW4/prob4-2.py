#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 21:49:57 2014

@author: eric
"""
import matplotlib.pyplot as plt
from ode import coupledRK4
from ode import coupledRK2

a = .08
b = .7
c = .8
i = .325

def f(x,y,t):
    return x-1/3*x**3-y+i 
    
def g(x,y,t):
    return a*(x+b-c*y)


x,y,t = coupledRK4(f,g,.1,1000,x0 = 1, y0 = 1)

plt.close()
plt.subplot(211)
plt.plot(x,y)
plt.subplot(212)
plt.plot(t,x,t,y)
plt.show()    
