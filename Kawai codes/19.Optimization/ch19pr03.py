#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 19.4.1                                                     *
%*     filename: ch19pr03.m                                               *
%*     program listing number: 19.3                                       *
%*                                                                        *
%*     This program fits a Gaussian disitrbution to a noisy data set      *
%*     using genetic algorithm.                                           *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  04/02/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# sample experimental data set
N=13
x=[-1.98,-1.48,-1.00,-0.50,-0.02,0.51,1.02,1.52,1.99,2.52,3.00,3.52,3.99]
y=[ 0.10, 0.11, 0.31, 0.65, 1.28,1.79,2.13,1.92,1.73,0.70,0.31,0.14,0.14]
s=[ 0.14, 0.11, 0.15, 0.18, 0.18,0.12,0.11,0.10,0.16,0.15,0.10,0.17,0.18]

# theoretical data space
K=101
X=np.linspace(x[0],x[-1],K)

# control parameter for genetic algorithm
NP=2**10           # population size
MP=np.int(NP/2)   # half of the population
max_count=200     # stopping condition
max_generation=10000   # maximum generation

# parameter space (genes)
a_max=5.0; a_min=0.0
b_max=5.0; b_min=-5.0
c_max=5.0; c_min=0.01

# initial parameters (genes)
a=np.random.rand(NP)*(a_max-a_min)+a_min
b=np.random.rand(NP)*(b_max-b_min)+b_min
c=np.random.rand(NP)*(c_max-c_min)+c_min

# allocate arrays
ip=np.zeros(MP)
F=np.zeros(NP)
G=np.zeros((N,NP))

# graphics setting
plt.close('all')
movie = True
if movie:
    plt.figure(figsize=(6,5))
    plt.axis([x[0], x[-1], -0.5, 2.5])
    
# reset the counters
found=False
count=0
generation=0
Fmin=np.finfo(np.float64()).max  # some large number 

# genetic evolution begins here
while not(found):
    
    generation=generation+1
    if generation > max_generation:  # too many generations
        print('Max generation reached.  Terminated.')
        break
    
    # eveluate the fitness
    for i in range(0,NP):
        for j in range(0,N):
            G[j,i]=a[i]*np.exp(-(x[j]-b[i])**2/c[i])
        F[i]=np.sum((G[:,i]-y[:])**2/s[:])/N

    # sort the population based on thier fitness
    IX=np.argsort(F)
    A=a[IX]; B=b[IX]; C=c[IX]
    
    if movie:
        #plot current best fitting
        plt.clf()
        Y=A[0]*np.exp(-(X-B[0])**2/C[0])
        plt.plot(X,Y,'-r',linewidth=2)
        plt.errorbar(x,y,yerr=s,fmt='ok')
        plt.pause(0.001)
    
    # check if converged
    if F[IX[0]]>=Fmin:
        count=count+1
        if count > max_count:  # no more evolution
            found=True
            break
    else:
        count=0;               # reset dewelling counter
        Fmin=F[IX[0]]          # new lowest fitness
        a0=A[0]; b0=B[0]; c0=C[0]  # current best genes
        print('New low found: generation={0:d}, fitness={1:.15f}, a={2:f}, b={3:f}, c={4:f}'.format(generation,Fmin,a0,b0,c0))

    # find mating pairs from the survived population
    ip=np.random.permutation(MP) 
     
    # generating oggsprings
    for k in range(0,MP,2):
        g=np.random.rand(6)
        i=ip[k]
        j=ip[k+1]
        A[k+MP]  = A[i]+g[0]*(A[j]-A[i]+0.01)
        A[k+MP+1]= A[j]+g[1]*(A[i]-A[j]+0.01)
        B[k+MP]  = B[i]+g[2]*(B[j]-B[i]+0.01)
        B[k+MP+1]= B[j]+g[3]*(B[i]-B[j]+0.01)
        C[k+MP]=   C[i]+g[4]*(C[j]-C[i]+0.01)
        C[k+MP+1]= C[j]+g[5]*(C[i]-C[j]+0.01)
 
    # mutation rate (every 10 generations, mutation burst happens)
    if np.mod(generation,10)==0:
        mutation=0.8
    else:
        mutation=0.1
        
    # mutation
    rm=np.random.rand(NP)
    for i in range(1,NP):
        if rm[i]<mutation:
            A[i]=np.random.rand(1)*(a_max-a_min)+a_min
            B[i]=np.random.rand(1)*(b_max-b_min)+b_min;
            C[i]=np.random.rand(1)*(c_max-c_min)+c_min;
    
    # store genes of new population
    a[:]=A[:]; b[:]=B[:]; c[:]=C[:]

print('Best fit:{0:.15f},  a={1:f},b={2:f}, c={3:f}'.format(Fmin, a0, b0, c0))

# plot the best fitting
if not(movie):
    plt.figure(figsize=(6,5))
    plt.axis([x[1], x[-1], -0.5, 2.5])
    
plt.clf()
Y=a0*np.exp(-(X-b0)**2/c0)
plt.plot(X,Y,'-r',linewidth=2)
plt.errorbar(x,y,yerr=s,fmt='ok')
plt.xlabel('x',fontsize=14)
plt.ylabel('f(x)',fontsize=14)
plt.show()
