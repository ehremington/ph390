#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 11.2                                                       *
%*     filename: ch11pr02.py                                              *
%*     program listing number: 11.2                                       *
%*                                                                        *
%*     This program calculates the fourier transform of the second order  *
%*     derivative.                                                        *
%*                                                                        *
%*     Uses Numpy FFT package                                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/18/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# control parameters
N=1024

# position domain
X=50.  # size of position domain
dx=X/N   # position interval
xmin=-X/2.
xmax=xmin+dx*(N-1)
x0=np.linspace(xmin,xmax,N)      # position space in natural order
x1=np.linspace(0.0,dx*(N-1),N)   # position space in fft order

# momentum domain
K=2.0*np.pi*N/X  # size of the momentum domain
dk=K/N           # resolution in momentum space
kmin=-K/2.
kmax=kmin+dk*(N-1)
k0=np.linspace(kmin,kmax,N)     # momentum space in natural order
k1=np.fft.fftshift(k0)          # momentum space in shifted order

f0=np.exp(-x0**2)/np.sqrt(np.pi) # input in natural order
f1=np.fft.fftshift(f0)           # input in fft order

# FFT from position to momentum domain
F1=np.fft.fft(f1)   # prefactor not need this time
F1=-F1*(k1**2)      # because we perform both forward
g1=np.fft.ifft(F1)  # and inverse transformation.
g0=np.fft.ifftshift(g1) # output in natural order

g2=np.exp(-x0**2)/np.sqrt(np.pi)*(4.0*x0**2-2.0)  # analytic FT

plt.figure(figsize=(6,5))
plt.plot(x0,np.real(g0),'-or',label='FFT')
plt.plot(x0,g2,'-b',label='exact')
plt.axis([-6., 6., -1.4, 0.8])
plt.xlabel('x',fontsize=14)
plt.ylabel('g(x)',fontsize=14)
plt.legend(loc=4)
plt.show()