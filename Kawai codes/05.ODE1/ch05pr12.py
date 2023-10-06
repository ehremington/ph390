#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 5.5.6                                                      *
%*     filename: ch05pr12.py                                              *
%*     program listing number: 5.12                                       *
%*                                                                        *
%*     This program calculate the trajectory of a particle scattered by   *
%*     Yukawa potential using he Valet method.                            *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt


E=np.float(input('Enter Energy = '))

# control parameter
bmax=3
N=10
db=bmax/N
h=0.01
Kmax=10001
b=np.linspace(-bmax,bmax,2*N+1)
theta=np.zeros(2*N+1)
x=np.zeros(Kmax)
y=np.zeros(Kmax)
vx=np.zeros(Kmax)
vy=np.zeros(Kmax)
plt.ioff()
plt.figure(figsize=(12,5))

plt.subplot(1,2,1);

for i in range(0,2*N+1):
    
    # initial conditions
    x[0]=-10; y[0]=b[i]; vx[0]=np.sqrt(2.0*E); vy[0]=0.0;
    
    # first Euler step
    r = np.sqrt(x[0]**2+y[0]**2)
    Fx=x[0]/r**2 * (1.0/r+1.0) * np.exp(-r)
    Fy=y[0]/r**2 * (1.0/r+1.0) * np.exp(-r)
    x[1]=x[0]+vx[0]*h+Fx*h**2/2.0
    y[1]=y[0]+vy[0]*h+Fy*h**2/2.0
    
   # Verlet method
    n=1
    while r <11.0:
        r = np.sqrt(x[n]**2+y[n]**2)
        Fx=x[n]/r**2 * (1.0/r+1.0) * np.exp(-r)
        Fy=y[n]/r**2 * (1.0/r+1.0) * np.exp(-r)
        x[n+1]=2.0*x[n]-x[n-1]+Fx*h**2
        y[n+1]=2.0*y[n]-y[n-1]+Fy*h**2
        n+=1
    
    # final velocity
    vfx=(x[n]-x[n-2])/(2.0*h)
    vfy=(y[n]-y[n-2])/(2.0*h)
    
    # scattering angle
    theta[i] = np.arccos((vx[0]*vfx+vy[0]*vfy)/(2.0*E))
    
    # plot the trajctory
    plt.plot(x[0:n],y[0:n],'-k')


# draw a target atom
plt.plot(0,0,'or')

plt.xlabel('x');
plt.ylabel('y');

# plot scattering angle
plt.subplot(1,2,2)
plt.plot(b,theta)
plt.xlabel('Impact Parameter')
plt.ylabel('Scattering Angle theta')
plt.show()


 