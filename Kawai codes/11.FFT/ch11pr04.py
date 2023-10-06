#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 11.4.4                                                     *
%*     filename: ch11pr04.py                                              *
%*     program listing number: 11.4                                       *
%*                                                                        *
%*     This program calculates the probability distribution of a harmonic *
%*     oscillator in the momentum spapce.                                 *
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
X=100.  # size of position domain
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
k1=np.linspace(0.0,dk*(N-1),N)          # momentum space in shifted order

# function in normal order
psix0=np.exp(-x0**2/2.0)/np.pi**(1./4.)
psix1=2.0*x0*np.exp(-x0**2/2.0)/np.pi**(1./4.)/np.sqrt(2.0)
psix2=(2.0*x0**2-1.0)*np.exp(-x0**2/2.0)/np.pi**(1./4.)/np.sqrt(2.0)

# Ground state (n=0)
f=np.fft.fftshift(psix0) # function in swapped order
g=np.fft.ifft(f)*X       # do not forget to multiply X
psik0=np.fft.ifftshift(g)
# 1st excited state (n=1)
f=np.fft.fftshift(psix1) # function in swapped order
g=np.fft.ifft(f)*X       # do not forget to multiply X
psik1=np.fft.ifftshift(g)
# 2nd excited state (n=2)
f=np.fft.fftshift(psix2) # function in swapped order
g=np.fft.ifft(f)*X       # do not forget to multiply X
psik2=np.fft.ifftshift(g)

N1=np.int(N/3)
N2=2*N1

plt.figure(figsize=(6,15))
plt.subplot(3,1,1)
plt.plot(k0[N1:N2],abs(psik0[N1:N2])**2)
plt.ylabel(r'$|\psi_0(k)|^2$')

plt.subplot(3,1,2)
plt.plot(k0[N1:N2],abs(psik1[N1:N2])**2)
plt.ylabel(r'$|\psi_1(k)|^2$')

plt.subplot(3,1,3)
plt.plot(k0[N1:N2],abs(psik2[N1:N2])**2)
plt.ylabel(r'$|\psi_2(k)|^2$')

plt.xlabel(r'$k$')
plt.show()


