#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 15.5.3                                                     *
%*     filename: ch15pr07.m                                               *
%*     program listing number: 15.7                                       *
%*                                                                        *
%*     This program simulatesa the surface growth using ballistic         *
%*     deposit model.                                                     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/26/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

N=20000  # number of particles to be deposited
L=200    #surface area

y=np.zeros(L,dtype=np.int)  # reset the height of the surface
x=np.random.randint(0,L,size=N)  # horizontal position of the particles

z=np.zeros(N)
w=np.zeros(N)
k=0
for i in range(0,N):
    y[x[i]]+=1  # ballistic growth
    
    # record the evolution of the growth after every 10 particles is deposited
    if np.mod(i,10)==0:
        z[k]=sum(y.astype(float))/L  # mean height
        w[k]=np.sqrt(sum((y.astype(float)-z[k])**2)/L)  # roughness
        k+=1

# calculate the height distribution
n1=min(y)  # lowest
n2=max(y)  # heighest
Ny=n2-n1+1
n=np.linspace(n1,n2,Ny)
h=np.zeros(Ny,dtype=np.int)

for i in range(0,L):
    j=y[i]-n1
    h[j]+=1
# normalization
h=h.astype(float)/sum(h)

# theoretical height distribution (Gaussian formula)
m=np.float(N)/L
s=np.float(N*(L-1))/L**2
g=np.exp(-(n-m)**2/(2.*s))/np.sqrt(2*np.pi*s)

# Figure 1: profile of the surface
plt.figure(figsize=(12,5))
X=np.linspace(1.0,L,L)
plt.bar(X,y,1.05,color='k')
plt.plot([0,L],[0,0],'-k',linewidth=4)  # draw the base line
plt.plot([0,L],[z[k-1],z[k-1]],'--r',label='Mean height')
plt.xlabel('x',fontsize=14)
plt.ylabel('height',fontsize=14)
plt.show()

# Figure 2: Evolution of the surface roughness
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(z[0:k],w[0:k],'-k',label='simulation')
plt.plot(z[0:k],np.sqrt(z[0:k]),'--r',label='theory')
plt.xlabel('height',fontsize=14)
plt.ylabel('surface roughness',fontsize=14)
plt.legend(loc=3)

# Figure 3: Heifht distribution
plt.subplot(1,2,2)
plt.bar(n,h,1.05,color='k')
plt.plot(n,g,'--r')

plt.xlabel('height',fontsize=14)
plt.ylabel('P(h)',fontsize=14)
plt.legend(loc=1)
plt.show()