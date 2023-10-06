# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 14:43:25 2014

@author: feynman
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
tf=30
h = 2**-12
#first derivatives
def dx(x,y,z):
    return s*(y-x)
def dy(x,y,z):
    return x*(p-z)-y
def dz(x,y,z):
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
    x_temp,y_temp,z_temp,t_temp = xn,yn,zn,tn
    xn = ode.RK4next(lambda x,t:dx(x,y_temp,z_temp),t_temp,x_temp,h)
    yn = ode.RK4next(lambda y,t:dy(x_temp,y,z_temp),t_temp,y_temp,h)  
    zn = ode.RK4next(lambda z,t:dz(x_temp,y_temp,z),t_temp,z_temp,h)
    tn = tn+h
    x_list.append(xn)
    y_list.append(yn)
    z_list.append(zn)  
    t_list.append(tn)


FFT = np.fft.fft(x_list)
FFT2 = FFT*np.conjugate(FFT)
freq = np.fft.fftfreq(t_list.__len__())   
t_array = np.asarray(t_list)
freq2 = 2*np.pi/t_array

#plt.close()
#plt.plot(freq,FFT2)
plt.plot(x_list,z_list)
plt.xlim((-5,5))
plt.yscale('log')
plt.show()   
