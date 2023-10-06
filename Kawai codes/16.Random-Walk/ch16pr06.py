#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#**************************************************************************
#*     Section 16.4.3                                                     *
#*     filename: ch16pr06.py                                              *
#*     program listing number: 16.6                                       *
#*                                                                        *
#*     This program simulates the Parrondo game.                          *
#*                                                                        *
#*     Programed by Ryoichi Kawai for Computational Physics Course.       *
#*     Last modification:  03/04/2017.                                    *
#**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# control parameters
e=0.005
pA=1./2.-e
qA=1.-pA
pB=3./4.-e
qB=1.-pB
pC=1./10.-e
qC=1.-pC
N=50000
M=100

# Game A
x=np.zeros(N)
y=np.zeros(M+1)
t=np.linspace(0,M,M+1)
for i in range(1,M+1):
    r=np.random.rand(N)
    x[r<pA]=x[r<pA]+1
    x[r>=pA]=x[r>=pA]-1
    y[i]=x.sum()/N

plt.figure(figsize=(6,5))
plt.plot(t,y,'-b',label='Game A alone',linewidth=2)

# Game B
x=np.zeros(N)
p=np.zeros(N)

for i in range(1,M+1):
    r=np.random.rand(N)
    p[:]=pB
    k=np.mod(x,3)==0
    p[k]=pC
    x[r<p]=x[r<p]+1
    x[r>=p]=x[r>=p]-1
    y[i]=x.sum()/N

plt.plot(t,y,'-g',label='Game B alone',linewidth=2)

# Alternating Game A and B
x=np.zeros(N)
p=np.zeros(N)

for i in range(1,M+1):
    r=np.random.rand(N)
    if np.mod(i,4)<2:
        p[:]=pA
    else:
        p[:]=pB
        k=np.mod(x,3)==0
        p[k]=pC

    x[r<p]=x[r<p]+1
    x[r>=p]=x[r>=p]-1
    y[i]=x.sum()/N

plt.plot(t,y,'-r',label='Game A & B',linewidth=2)

plt.legend(loc=2)
plt.xlabel('# of plays',fontsize=14)
plt.ylabel('Gain',fontsize=14)
plt.show()


