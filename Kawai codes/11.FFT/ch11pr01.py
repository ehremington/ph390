#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 11.1                                                       *
%*     filename: ch11pr01.py                                              *
%*     program listing number: 11.1                                       *
%*                                                                        *
%*     This program calculates the fourier transform of Gaussian.         *
%*                                                                        *
%*     Uses Numpy fft package                                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/18/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# control parameters
N=1024

# Setting time domain
T=50.   # size of time domain
dt=T/N  # time interval
tmin=-T/2.
tmax=tmin+dt*(N-1)
t0=np.linspace(tmin,tmax,N)    # original time
t1=np.linspace(0.0,dt*(N-1),N) # fft ordered time

# setting frequency domain
W=2.0*np.pi*N/T # size of the frequency domain
dw=W/N          # frequency resolution
wmin=-W/2.0
wmax=wmin+dw*(N-1)
w0=np.linspace(wmin,wmax,N)     # frequency domain we want
w1=np.linspace(0.0,dw*(N-1),N)  # fft frequency domain

# original Gaussian function in time domain
f0=np.exp(-t0**2)/np.sqrt(np.pi)
# function in fft order
f1=np.fft.fftshift(f0)

# FFT from time to frequency domain in  ffto order
g1=np.fft.ifft(f1)*T  # do not forget to multiply T

# Get the normal frequenx=cy domain order
g0=np.fft.ifftshift(g1)

# analytic FT
g2=np.exp((-w0**2)/4.0)

plt.figure(figsize=(12,10))
plt.subplot(2,2,1)
plt.plot(t0,f0,'-r',label=r'original input',linewidth=2)
plt.plot([-T/2,-T/2],[0,1],'--k')
plt.plot([T/2,T/2],[0,1],'--')
plt.axis([-T/2*1.05, T*1.05, 0.0, 1.0])
plt.legend(loc=2)
plt.xlabel(r'$t$',fontsize=14)
plt.ylabel(r'$f(t)$',fontsize=14)

plt.subplot(2,2,3)
plt.plot(t1,f1,'-r',label='fft input',linewidth=2)
plt.plot([0,0],[0,1],'--k')
plt.plot([T,T],[0,1],'--k')
plt.axis([-T/2*1.05, T*1.05, 0.0, 1.0])
plt.legend(loc=2)
plt.xlabel(r'$t$',fontsize=14)
plt.ylabel(r'$f(t)$',fontsize=14)

plt.subplot(2,2,2)
plt.plot(w1,np.real(g1),'-r',label='fft output',linewidth=2)
plt.plot([0,0],[0,1.3],'--k')
plt.plot([W,W],[0,1.3],'--k')
plt.axis([-W/2*1.05, W*1.05, 0.0, 1.3])
plt.legend(loc=2)
plt.xlabel(r'$\omega$',fontsize=14)
plt.ylabel(r'$\tilde{f}(\omega)$',fontsize=14)

plt.subplot(2,2,4)
plt.plot(w0,g0,'-r',label='desired output',linewidth=2)
plt.plot([-W/2,-W/2],[0,1.3],'--k')
plt.plot([W/2,W/2],[0,1.3],'--k')
plt.axis([-W/2*1.05, W*1.05, 0.0, 1.3])
plt.legend(loc=2)
plt.xlabel(r'$\omega$',fontsize=14)
plt.ylabel(r'$\tilde{f}(\omega)$',fontsize=14)
plt.show()

plt.figure(figsize=(6,5))
plt.plot(w0,g0,'-or',label='FFT',linewidth=2)
plt.plot(w0,g2,'-k',label='Exact',linewidth=2)
plt.xlabel(r'$\omega$',fontsize=14)
plt.ylabel(r'$\tilde{f}(\omega)$',fontsize=14)
plt.show()
