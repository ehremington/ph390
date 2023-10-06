#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 15.5.2                                                     *
%*     filename: ch15pr06.py                                              *
%*     program listing number: 16.6                                       *
%*                                                                        *
%*     This program generates distribution of particles under gravity     *
%*     thermal diffusion.
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/26/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

N=1000  # number of particles

 # horizontal position = uniform
x=np.random.rand(N)
# vertical position = exponential
y=np.random.rand(N) # vertical position = exponential
z=-np.log(y)
zmax=np.max(z)+1.


plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot([-0.05,-0.05],[-0.05, zmax],'-k',linewidth=2)
plt.plot([-0.05, 1.05],[-0.05,-0.05],'-k',linewidth=2)
plt.plot([ 1.05, 1.05],[-0.05, zmax],'-k',linewidth=2)
plt.plot(x,z,'.');
plt.axis([-0.1, +1.1, -0.5, zmax])
plt.xlabel('x',fontsize=14)
plt.ylabel('z',fontsize=14)
plt.subplot(1,2,2)
plt.hist(z,2*np.int(zmax),normed=1)
Z=np.linspace(0.0,zmax,201)
P=np.exp(-Z)
plt.plot(Z,P,'--r',linewidth=2)
plt.xlabel('z',fontsize=14)
plt.ylabel('probability density',fontsize=14)
plt.show()