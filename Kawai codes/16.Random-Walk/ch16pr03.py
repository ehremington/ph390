#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Exercise 17.3                                                      *
%*     filename: ch17pr03.m                                               *
%*     program listing number: 17.3                                       *
%*                                                                        *
%*     This program simulates the two-dimensional random walk.            *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/15/2014.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

M=100000  # number of particles
N=1000    # number of steps

x=np.zeros(M)  # initial position
y=np.zeros(M)

u=np.zeros((N+1,2))
w=np.zeros((N+1,2))
s=np.zeros(N+1)

for i in range(1,N+1):

    g=np.random.random_integers(1,4,M)  # pick one of four directions to jump
    
    # jump
    # The following expression is easy to write but slow in execution
    x[g==1]+=1
    x[g==2]+=-1
    y[g==3]+=1
    y[g==4]+=-1
    
    # record two sample trajectories
    u[i,0]=x[0]
    w[i,0]=y[0]
    u[i,1]=x[1]
    w[i,1]=y[1]

    # mean square displacenent
    s[i]=sum(x**2+y**2)/M

fig1, ax=plt.subplots(figsize=(6,6))
plt.plot(u[:,0],w[:,0],'-b')
plt.plot(u[:,1],w[:,1],'-g')
R=np.sqrt(s[N])
c=plt.Circle((0, 0), R, color='r',fill=False)
ax.add_artist(c)
ax.axis('equal')
L=R*1.5
ax.axis([-L, L, -L, L])
plt.xlabel('x',fontsize=14)
plt.ylabel('y',fontsize=14)
plt.legend(loc=1)

fig2, bx=plt.subplots(figsize=(6,6))

plt.plot(x[0:2000],y[0:2000],'.')
c=plt.Circle((0, 0), R, color='r',fill=False,linewidth=2)
bx.add_artist(c)
bx.axis('equal') 
bx.axis([-100, 100, -100, 100]) 

plt.xlabel('x',fontsize=14)
plt.ylabel('y',fontsize=14)

plt.show()