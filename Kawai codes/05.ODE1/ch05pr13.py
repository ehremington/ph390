#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 5.5.7                                                      *
%*     filename: ch05pr13.m                                               *
%*     program listing number: 5.13                                       *
%*                                                                        *
%*     This program calculate the trajectory of a double pendulum by      *
%*     the 4th order Runge-Kutta method.                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# system parameters 
m1=2; m2=1; L1=1; L2=2; g=9.8;

def dp(q1,q2,w1,w2):
    global m1, m2, L1, L2, g
    A = (m1+m2)*L1
    B = m2*L2*np.cos(q1-q2)
    C = m2*L1*np.cos(q1-q2)
    D = m2*L2
    E = -m2*L2*w2**2*np.sin(q1-q2)-(m1+m2)*g*np.sin(q1)
    F = m2*L1*w1**2*np.sin(q1-q2)-m2*g*np.sin(q2)
    Y1=(D*E-B*F)/(A*D-B*C)
    Y2=(A*F-C*E)/(A*D-B*C)
    return [Y1, Y2]

# control parameter
tmax=100.0
N=10000
h=tmax/N

# allocate araay
t=np.linspace(0,tmax,N+1)
q1=np.zeros(N+1)
q2=np.zeros(N+1)
w1=np.zeros(N+1)
w2=np.zeros(N+1)

# initial conditions
q1[0]=1.5
q2[0]=3.0
w1[0]=0.0
w2[0]=0.0

V = -(m1+m2)*g*L1*np.cos(q1[0])-m2*g*L2*np.cos(q2[0])
T = m1*L1**2*w1[0]**2/2.0+m2*(L1**2*w1[0]**2+L2**2*w2[0]**2\
    + 2*L1*L2*w1[0]*w2[0]*np.cos(q1[0]-q2[0]))/2.0
E = T+V
print("Initial Energy ={0:16.8e}".format(E))

for n in range(0,N):
    
    # 4th-order Runge-Kutta
    dotw=dp(q1[n],q2[n],w1[n],w2[n])
    kw11=dotw[0]
    kw21=dotw[1]
    kq11=w1[n]
    kq21=w2[n]
    
    w1m = w1[n]+kw11*h/2.0
    w2m = w2[n]+kw21*h/2.0
    q1m = q1[n]+kq11*h/2.0
    q2m = q2[n]+kq21*h/2.0
    
    dotw=dp(q1m,q2m,w1m,w2m)
    kw12 = dotw[0]
    kw22 = dotw[1]
    kq12 = w1m
    kq22 = w2m
    
    w1m = w1[n]+kw12*h/2.0
    w2m = w2[n]+kw22*h/2.0
    q1m = q1[n]+kq12*h/2.0
    q2m = q2[n]+kq22*h/2.0

    dotw=dp(q1m,q2m,w1m,w2m)
    kw13 = dotw[0]
    kw23 = dotw[1]
    kq13 = w1m
    kq23 = w2m
    
    w1f = w1[n]+kw13*h
    w2f = w2[n]+kw23*h
    q1f = q1[n]+kq13*h
    q2f = q2[n]+kq23*h
    
    dotw=dp(q1f,q2f,w1f,w2f)
    kw14 = dotw[0]
    kw24 = dotw[1]
    kq14 = w1f
    kq24 = w2f
    
    q1[n+1]=q1[n]+(kq11+2.0*(kq12+kq13)+kq14)*h/6.0
    q2[n+1]=q2[n]+(kq21+2.0*(kq22+kq23)+kq24)*h/6.0
    w1[n+1]=w1[n]+(kw11+2.0*(kw12+kw13)+kw14)*h/6.0
    w2[n+1]=w2[n]+(kw21+2.0*(kw22+kw23)+kw24)*h/6.0

V = -(m1+m2)*g*L1*np.cos(q1[n+1])-m2*g*L2*np.cos(q2[n+1]);
T = m1*L1**2*w1[n+1]**2/2.0+m2*(L1**2*w1[n+1]**2+L2**2*w2[n+1]**2
        +2.0*L1*L2*w1[n+1]*w2[n+1]*np.cos(q1[n+1]-q2[n+1]))/2.0
E=T+V
print("  Final Energy ={0:16.8e}".format(E))

# plot angular coordinates
plt.ioff()
plt.figure(figsize=(12,5))
plt.axes().set_aspect('equal')
plt.subplot(1,2,1)
plt.plot(t,q1,'-b',label=r'$\theta_1$')
plt.plot(t,q2,'-r',label=r"$\theta_2$")
plt.xlabel("t")
plt.ylabel(r"$\theta$")
plt.legend(loc=3)

# trajectory of the second bob in xy coordiates
plt.subplot(1,2,2)
x2=L1*np.sin(q1)+L2*np.sin(q2)
y2=-L1*np.cos(q1)-L2*np.cos(q2)
plt.plot(x2,y2,'-k')
plt.xlabel('x')
plt.ylabel('y')

plt.show()