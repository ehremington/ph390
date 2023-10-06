#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section  9.3.1                                                     *
%*     filename: ch09pr02.py                                               *
%*     program listing number: 9.2-1                                      *
%*                                                                        *
%*     This program finds a fixed point of Maxwell-Bloch equation         *
%*     by newton-raphson method and Gaussian elimination.                 *
%*                                                                        *
%*     Use function gauss.m                                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/15/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt


# system parameters
g1=0.1; g2=2.0; g3=3.0
k1=0.25; k2=0.2; k3=1.0
lam=5.0

def Jacob(x):
    J = [[-g1,k1,0.0],[k2*x[2],-g2,k2*x[1]],[-k3*x[1],-k3*x[0],-g3]]
    return J
    
def Func(x):
    f= [-g1*x[0]+k1*x[1],-g2*x[1]+k2*x[0]*x[2],-g3*(x[2]-lam)-k3*x[0]*x[1]]
    return f
            
# control parameters
alpha = 0.5
tol = 1e-4

nmax=1000
y=np.zeros((3,nmax+1))
J=np.zeros((3,3))
f=np.zeros(3)
B=np.zeros(3)

#initial guess
x=np.array([2.0,2.0,5.0])
y[:,0]=x
#y[:,0]=[x.item(0),x.item(1),x.item(2)]

n=0
err=tol+1.0
J=np.array(Jacob(x))
f=np.array(Func(x))

while err > tol and n<nmax:
    B=np.linalg.solve(J,f)
    x=y[:,n] - alpha*B
    y[:,n+1]=[x.item(0),x.item(1),x.item(2)]
    J=np.array(Jacob(x))
    f=np.array(Func(x))
    err = np.dot(f,f)
    n+=1

print('E={0:f},  P={1:f},  D={2:f}'.format(x.item(0),x.item(1),x.item(2)))
print('err=',err)

t=np.linspace(0,n,n+1)
plt.figure()
plt.plot(t,y[0,0:n+1],'-ok',label='E')
plt.plot(t,y[1,0:n+1],'-ob',label='P')
plt.plot(t,y[2,0:n+1],'-or',label='D')
plt.xlabel('iteration')
plt.ylabel('lasing state')
plt.xlim(0.,9.)
plt.ylim(0.,6.)
plt.legend(loc=1)
plt.show()
