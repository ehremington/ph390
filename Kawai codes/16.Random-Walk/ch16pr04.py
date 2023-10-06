#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#**************************************************************************
#*     Section 16.4.1                                                     *
#*     filename: ch16pr04.py                                              *
#*     program listing number: 16.4                                       *
#*                                                                        *
#*     This program simulates the two-dimensional diffusion limited       *
#*     aggregates.                                                        *
#*                                                                        *
#*     Programed by Ryoichi Kawai for Computational Physics Course.       *
#*     Last modification:  03/04/2017.                                    *
#**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

anim=True

# number of particles
N=10000

L=601 # size of the square space
L0=np.int(L/2) # center of the square

# set array size
x=np.zeros(N,dtype=np.int)
y=np.zeros(N,dtype=np.int)
A=np.zeros((L,L),dtype=np.int)

R=5.         # inner circle
R_max=3.*R   # outer circle

# seed particle
A[L0,L0]=1
x[0]=L0
y[0]=L0

plt.ion()
fig, ax = plt.subplots(figsize=(8,8))
ax.set_xlim([L0-100,L0+100])
ax.set_ylim([L0-100,L0+100])

c=plt.Circle((x[0],y[0]), 0.5, color='b')
ax.add_artist(c)


for n in range(1,N):
    # random point on the inner circle
    theta=np.random.rand(1)*2.0*np.pi
    x[n]=np.int(R*np.cos(theta))+L0
    y[n]=np.int(R*np.sin(theta))+L0
    
    # diffusion
    found=False 
    while not(found):
        p=np.random.rand(1)
        if p<1./4.:
            x[n]+=1
        elif p<1./2.:
            x[n]+=-1
        elif p<3./4.:
            y[n]+=1
        else:
            y[n]+=-1

        r=np.sqrt(np.float((x[n]-L0)**2+(y[n]-L0)**2))
        if r>R_max:  # out of bound - restart
            theta=np.random.rand(1)*2.0*np.pi
            x[n]=np.int(R*np.cos(theta))+L0
            y[n]=np.int(R*np.sin(theta))+L0
        elif r<R:
            # hit the cluster?
            if A[x[n]+1,y[n]]+A[x[n]-1,y[n]]+A[x[n],y[n]+1]+A[x[n],y[n]-1]>0:
                found=True
                A[x[n],y[n]]=1
                c=plt.Circle((x[n],y[n]), 0.5, color='b',fill=False)
                ax.add_artist(c)
                if anim:
                    plt.pause(0.0001)
                if R<r+5:
                    R=r+5. # adjust inner circle radius
                R_max = 3*R     # adjust outer circle radius
                if R>L0:
                    plt.xlabel('x',fontsize=14)
                    plt.ylabel('y',fontsize=14)
                    plt.show()
                    exit('Out of Range')


