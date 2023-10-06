#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 12.3.2                                                     *
%*     filename: ch12pr08.m                                               *
%*     program listing number: 12.8                                       *
%*                                                                        *
%*     This program finds the peak position and life-time broadening      *
%*     of atomic emmision spectrum.                                       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# Generate a sample noisy data
x=[-2.01,-1.47,-0.97,-0.52,-0.04,0.52,0.99,1.53,2.03,2.51,2.96,3.47,4.02]
y=[0.28,0.57,0.62,0.68,1.26,1.29,1.57,1.11,0.91,0.94,0.65,0.80,0.31]
s=[0.10,0.11,0.17,0.06,0.15,0.11,0.15,0.10,0.11,0.14,0.16,0.18,0.15]

x=np.array(x)
y=np.array(y)
s=np.array(s)
N=x.size

# control parametrs
alpha=1.e-2


# initial guess
K=2000
M=3
lam=np.zeros((K,M))
chi2=np.zeros(K)

lam[0,:]=np.array([1.,0.,1.])

# Gauss-Newton iteration
n=0
# evaluate initial chi sqaure
b=np.zeros(N)
J=np.zeros((N,M))

for i in range(0,N):
   F=lam[n,0]/((x[i]-lam[n,1])**2+lam[n,2])
   b[i]=y[i]-F

b=b/s
chi2[n]=(b**2).sum()


found=False
while not(found):
    # construct Jacobian and vectors.
    for i in range(0,N): 
        F=lam[n,0]/((x[i]-lam[n,1])**2+lam[n,2])
        J[i,0]=F/lam[n,0]
        J[i,1]=2.*lam[n,0]*(x[i]-lam[n,1])\
              /((x[i]-lam[n,1])**2+lam[n,2])**2
        J[i,2]=-lam[n,0]/((x[i]-lam[n,1])**2+lam[n,2])**2
        b[i]=y[i]-F

    # Take into account error bar
    b=b/s
    for i in range(0,M):
       J[:,i]=J[:,i]/s
   
    # Solve the equation
    n+=1
    A=np.dot(J.transpose(),J)
    c=np.dot(J.transpose(),b)
    d=np.linalg.solve(A,c)
    
    # update parameter values
    lam[n,:]=lam[n-1,:]+alpha*d
    # evaluate chi sqaure
    chi2[n]=(b**2).sum()    

    # if chi square goes up stop
    if chi2[n]-chi2[n-1] > 0:
        found=True

        
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
z=np.linspace(-4.0,6.0,101)
L=z.size
w=np.zeros(L)

for i in range(0,L):
    w[i]=lam[n,0]/((z[i]-lam[n,1])**2+lam[n,2])

plt.plot(z,w)
plt.errorbar(x,y,yerr=s,fmt='o')
plt.xlabel('x',fontsize=14)
plt.ylabel('f(x)',fontsize=14)

plt.subplot(1,2,2)
plt.semilogy(np.linspace(1,n,n),chi2[0:n])
plt.xlabel('iteration',fontsize=14)
plt.ylabel(r'$\chi^2$',fontsize=14)
plt.show()
