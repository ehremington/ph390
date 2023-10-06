#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 15.5.3                                                     *
%*     filename: ch15pr09.m                                               *
%*     program listing number: 15.9                                       *
%*                                                                        *
%*     This program simulatesa the surface growth using ballistic         *
%*     deposit model with overhangs.                                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/26/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

movie=False  # Set this to True to show real time growth (very slow)

N=20000
L=200
fig, ax=plt.subplots()

ax.set_xlim([0,L])
ax.set_ylim([0,N/L*3])


y=np.zeros(L,dtype=np.int)  # reset the height of the surface
x=np.random.randint(0,L,size=N)  # horizontal position of the particles

z=np.zeros(N)
w=np.zeros(N)
k=0
for i in range(0,N):
    
    # lateral diffusion
    j0=x[i]
    j1=np.mod(j0-1,L) # left neighbor
    j2=np.mod(j0+1,L) # right neighbor
    if y[j0]<y[j1] or y[j0]<y[j2]:
        y[j0]=np.max( (y[j1],y[j2]) )  # stick to the next site
    else:
        y[j0]+=1                 # regular deposition
        
    # draw the particle
    c=plt.Circle((j0,y[j0]), 0.5, color='b')
    ax.add_artist(c)
    if movie:
        plt.pause(0.0001)
    
    # record the evolution of the growth after every 10 particles is
    # deposited    
    if np.mod(i,10)==0:
        z[k]=sum(y.astype(float))/L  # mean height
        w[k]=np.sqrt(sum((y.astype(float)-z[k])**2)/L)  # roughness
        k+=1
        
ax.plot([0,L],[z[k-1],z[k-1]],'--r')
plt.show()
