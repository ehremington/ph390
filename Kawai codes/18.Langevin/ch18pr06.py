#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 18.2.2                                                     *
%*     filename: ch18pr06.m                                               *
%*     program listing number: 18.6                                       *
%*                                                                        *
%*     This program simulates a stochastic resonace.                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/26/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt


# system parameters
# To see stochastic resonance, use the parameter valesu
# A=1.0; D=1.2; U0=10; w=0.5;
#
A=1.0
D=1.2
U0=10
w=0.5

def U(x,t):
    return U0*(x**4/4.0-x**2/2.0)+0.05-A*np.cos(w*t)*x

def F(x,t):
    return U0*(-x**3+x)+A*np.cos(w*t)

# constrol parameters
movie=False
tau=100.0
N=2**14
dt=tau/N
ds=np.sqrt(2.0*D*dt)

t=np.linspace(0.0,dt*(N-1),N)
x=np.zeros(N)

L=101
q=np.linspace(-2.0,2.0,L)
U1=U0*(q**4/4.0-q**2/2.0)

plt.close('all')
if movie:
    fig1=plt.figure(figsize=(4,7))

    plt.plot(q,U1-A*np.cos(w*t[0])*q,'-b')
    y=U(x[0],t[0])
    circle=plt.Circle((x[0],y),0.05,fc='k')
    plt.gca().add_patch(circle)
    plt.axis([-2, 2, -U0/2, 4])
#    plt.axis('equal')
    plt.pause(0.001)

g=np.random.randn(N)

for i in range(0,N-1):

    f1=F(x[i],t[i])
    x[i+1]=x[i] + f1*dt + g[i]*ds
    f2=F(x[i+1],t[i+1])
    x[i+1]=x[i] + (f1+f2)*dt/2.0 + g[i]*ds
    if movie :
        plt.clf()

        plt.plot(q,U1-A*np.cos(w*t[i+1])*q,'-b')

        plt.axis([-2, 2, -U0/2, 4])
#        plt.axis('equal')
        y=U(x[i+1],t[i+1])
        circle=plt.Circle((x[i+1],y),0.1,fc='k')
        plt.gca().add_patch(circle)
        plt.pause(0.001)

        
plt.figure(figsize=(6,5))
plt.plot(t,x)
plt.plot(t,np.cos(w*t),'--r')
plt.xlabel('t',fontsize=14)
plt.ylabel('x(t)',fontsize=14)
plt.show()

plt.figure(figsize=(6,5))
z=np.fft.fft(x)*tau/N
u=np.linspace(0.0,2.0*np.pi/tau*(N-1),N)
plt.plot(u,abs(z)**2)
plt.xlim([0,3])
plt.ylabel('Power spectrum',fontsize=14)
plt.xlabel(r'$\omega$',fontsize=14)
plt.plot([0.5,0.5],[0,900],'--k');
plt.show()

