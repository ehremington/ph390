#!/usr/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 16:31:16 2014

@author: eric
"""


import ph762.ode as ode
import numpy as np
import matplotlib.pyplot as plt
#parameter values
s=10
b=8/3
p=28

#initial conditions
x0=10
y0=10
z0=10
t0=0
tf=1024
h = 2**-7
#first derivatives
def dx(x,y,z,t):
    return s*(y-x)
def dy(x,y,z,t):
    return x*(p-z)-y
def dz(x,y,z,t):
    return x*y-b*z
    
x_list = [x0]
y_list = [y0]
z_list = [z0]
t_list = [t0]
xn=x0
yn=y0
zn=z0
tn=t0

while tn<tf:
    xn,yn,zn = ode.cRK4next(dx,dy,dz,zn,yn,xn,tn,h)
    tn = tn + h
    x_list.append(xn)
    y_list.append(yn)
    z_list.append(zn)
    t_list.append(tn)



FFT = np.fft.fft(x_list)
FFT2 = np.abs(FFT)**2
freq = np.fft.fftfreq(t_list.__len__(),h)   
plt.close()
plt.plot(freq,FFT2,'.')
#plt.plot(x_list,z_list)
plt.xlim((0,20))
plt.yscale('log')
plt.show()   