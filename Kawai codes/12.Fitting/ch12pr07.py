#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 12.3.1                                                     *
%*     filename: ch12pr07.py                                               *
%*     program listing number: 12.7                                       *
%*                                                                        *
%*     This program finds the activation energy of a reaction from a data *
%*     set using the linear regression.                                   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                     *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

k=8.617e-5; # Boltzmann constant [eV/K]

T=[200., 220., 240., 260., 280., 300., 320., 340., 360., 380., 400.]
f=[0.471,0.515,0.576,0.639,0.734,0.742,0.833,0.830,0.932,0.918,0.939]
f=np.array(f)
T=np.array(T)
N=f.size
x=1./(k*T[::-1])   
y=np.log(f[::-1])

# Linear regression
F=y.sum()
X=x.sum()
X2=(x**2).sum()
XF=(x*y).sum()
b=(F*X2-X*XF)/(N*X2-X**2)
a=(N*XF-X*F)/(N*X2-X**2)

g=a*x+b
A=np.exp(b)
z = A*np.exp(a/(k*T))

print('\nActivation Energy = {0:8.4f} eV\n'.format(np.abs(a)))

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(x,g,'-r')
plt.plot(x,y,'ok')
plt.xlabel(r'$\beta$',fontsize=14)
plt.ylabel(r'$\log\, k$',fontsize=14)

plt.subplot(1,2,2)
plt.plot(T,f,'ok')
plt.plot(T,z,'-r')
plt.xlabel('T',fontsize=14)
plt.ylabel('k',fontsize=14)
plt.show()