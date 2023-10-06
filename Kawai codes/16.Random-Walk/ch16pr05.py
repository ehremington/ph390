#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#**************************************************************************
#*     Section 17.4.2                                                     *
#*     filename: ch17pr05.m                                               *
#*     program listing number: 17.5                                       *
#*                                                                        *
#*     This program simulates the growth of dendrite.                     *
#*                                                                        *
#*     Programed by Ryoichi Kawai for Computational Physics Course.       *
#*     Last modification:  02/15/2014.                                    *
#**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

anim=True

# size of the systrem
N=5000 
Lx=100 
Ly=100 

plt.ion()
fig, ax = plt.subplots(figsize=(6,6))
plt.axis('equal')
ax.set_xlim([0,Lx])
ax.set_ylim([0,Ly])

# initial setting
x=np.random.random_integers(0,Lx-1,N)  #  horizontal position (uniform random)
y=np.zeros(N,dtype=np.int) 
A=np.zeros((Lx,2*Ly),dtype=np.int) 
A[:,0]=1 
H=5 
ymax=5 

# bias in y direction
e=0.01 

# initial plots
plt.plot([0,Lx-1],[0.5,0.5],'-k',linewidth=2) 

# deposition process
for n in range(0,N):
    y[n]=H   # diffusion starts from here 
    
    found=False 
    
    # diffues in 2D space until it sticks to another.
    while not(found):
        p=np.random.rand(1) 
        if p<1./4.:
            x[n]=np.mod(x[n]+1,Lx)   # periodic boundary
        elif p<1./2.:
            x[n]=np.mod(x[n]-1,Lx)   # periodic boundary
        elif p<3./4.-e:  
            y[n]=y[n]+1 
            if y[n]>3*H:  # if it went to high, start over.
               y[n]=H 
        else:
            y[n]=y[n]-1 

        if y[n]<H:
            i1=np.mod(x[n]-1,Lx) 
            i2=np.mod(x[n]+1,Lx) 
            if A[i1,y[n]]+A[i2,y[n]]+A[x[n],y[n]+1]+A[x[n],y[n]-1]>0:
                found=True 
                A[x[n],y[n]]=1 
                c=plt.Circle((x[n],y[n]), 0.5, color='b',fill=False)
                ax.add_artist(c)
                if anim:
                    plt.pause(0.0001)

                ymax=np.max([y[n],ymax])  # adjust the stating height
                H=5+ymax 
                if ymax>Ly-1:
                    plt.xlabel('x',fontsize=14)
                    plt.ylabel('y',fontsize=14)
                    plt.show()
                    exit('Out of Range')

