#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Exercise 17.2                                                      *
%*     filename: ch17pr02.py                                              *
%*     program listing number: 17.2                                       *
%*                                                                        *
%*     This program simulates two-dimensional Ising model using           *
%*     the Metropolis algorithm.                                          *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/17/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# ease all previous figures
plt.close('all')

# control parameters
L=32          # number of spins in one direction
LL=L*L        # total number of spins
N0=20000      # thermalization steps
N=100000      # number of Metropolis steps.
NS=N/1000     # number of samples.

# Define arrays
kmax=np.int(N/NS)
m=np.zeros(kmax)
sigma2=np.zeros(kmax)
E=np.zeros(kmax)

# Show animation if True
movie=False

# temperature
T=2.0        
beta = 1.0/T

# Initial configuration
if T>2.0:
    # Random (for high temperature)
    s=np.random.choice([1,-1],[L,L])
else:
    # Uniform (for low temperature)
    s=np.ones((L,L),dtype=np.int)

#  Show initial configuration
if movie:
    plt.figure(figsize=(6,6))
    plt.axis('equal')
    plt.axes(xlim=(-1, L), ylim=(-1, L))
    for j in range(0,L):
        for i in range(0,L):
            if s[i,j]==1:
                color='b'
            else:
                color='y'
            circle=plt.Circle((i,j),0.5,fc=color)
            plt.gca().add_patch(circle)
    plt.pause(0.0001)

# Begin Metropolis simulation    
k=0
for n in range(0,N+N0):
    # pick a site at random 
    i=np.random.randint(0,L)
    j=np.random.randint(0,L)

    # Evaluation of energy change
    i1=np.mod(i+1,L)
    i2=np.mod(i-1,L)
    j1=np.mod(j+1,L)
    j2=np.mod(j-1,L)
    ss = s[i1,j]+s[i2,j]+s[i,j1]+s[i,j2]
    dE = 2*ss*s[i,j]

# Flip spin based on Metropolis algorithm
    if np.exp(-beta*dE)> np.random.rand(1):
        s[i,j]=-s[i,j]

        # Show new configuration
        if movie:
            if s[i,j]==1:
                color='b'
            else:
                color='y'
            circle=plt.Circle((i,j),0.5,fc=color)
            plt.gca().add_patch(circle)
            plt.pause(0.0001)
            
    # Evaluate statistical quantities
    if n>N0 and np.mod(n,NS)==0:
        
        # mean and variance
        m[k] = np.real(s.sum())/LL
        sigma2[k] = (s**2).sum()/LL-m[k]**2 

        # total nergy
        h=0
        for j in range(0,L-1):
            for i in range(0,L-1):
                 h=h+s[i,j]*(s[i+1,j]+s[i,j+1])
        for i in range(0,L):
            h=h+s[i,L-1]*s[i,0]+s[L-1,i]*s[0,i]
        E[k]=-h

        k+=1

# plot magnetization
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
t=np.linspace(0,k-1,k)*NS
plt.plot(t,m[0:k],'-b',label='m(t)')

mu=sum(m[1:k])/k
plt.plot([0, N],[mu, mu],'--r',label='mean')
plt.xlim([0,N])
plt.ylim([-1.1,1.1])
plt.legend(loc=4)
plt.xlabel('steps')
plt.ylabel('magnetization')

# plot energy
plt.subplot(1,2,2)
Eavg=sum(E[0:k])/k
C=(sum(E[0:k]**2)/k-Eavg**2)/T**2/LL
plt.plot(t,E[0:k]/LL,label='energy');
plt.plot([0,N],[Eavg/LL, Eavg/LL],'--r',label='mean')
plt.xlim([0,N])
plt.ylim([-4,0])
plt.xlabel('steps')
plt.ylabel('Energy/spin')
plt.show()

# Show the final configulation
if not(movie):
    plt.figure(figsize=(6,6))
    plt.axis('equal')
    plt.axes(xlim=(-1, L), ylim=(-1, L))
    for j in range(0,L):
        for i in range(0,L):
            if s[i,j]==1:
                color='b'
            else:
                color='y'
            circle=plt.Circle((i,j),0.5,fc=color)
            plt.gca().add_patch(circle)
    plt.show()
    
# statistics
print('<m>={0:8.5f}'.format(mu))
print('<E>={0:8.5f}'.format(Eavg/LL))
print('<C>={0:8.5f}'.format(C))


