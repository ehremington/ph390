#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  8.11                                                      *
%*     filename: ch08pr10.py                                              *
%*     program listing number: 8.10                                       *
%*                                                                        *
%*     This program solves a 2x2 linear equation by the conjugate         *
%*     gradient minimization.                                             *
%*                                                                        *
%*     We use numpy matrix instead of array and column vectors as         *
%*     transpose of row vector.  Numpy treats it as (2x1) matrix.         *
%*     Because of that dotproduct becomes (1x1) matrix.  Use asscalar()   *
%*     in Numpy to correct it.
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

A=np.matrix([[4.,1.],[1.,3.]])
b=np.matrix([1.,2.]).transpose()
p=np.matrix([0.,0.]).transpose()
r=np.matrix([0.,0.]).transpose()
r1=np.matrix([0.,0.]).transpose()

tol = 1e-8;

# conjugate gradient method
kmax=1000
tol = 1e-8
x=np.matrix([1,0.5]).transpose() # starting point
u=np.zeros(kmax)
v=np.zeros(kmax)
u[0]=np.asscalar(x[0])
v[0]=np.asscalar(x[1])

r=b-A*x
p[:]=r[:]
rr=np.asscalar(r.transpose()*r)
n=0
while abs(rr)>tol and n<kmax:
    n+=1
    alpha = rr/np.asscalar(r.transpose()*A*r)
    x = x + alpha*p
    u[n]=np.asscalar(x[0])  # In numpy, a column vector must be 
    v[n]=np.asscalar(x[1])  # treated as matrix of (Nx1).
    r1=b-A*x
    rr1=np.asscalar(r1.transpose()*r1)
    beta=rr1/rr
    r[:]=r1[:]
    rr=rr1
    p[:]=r[:]+beta*p[:]


print('\nSolution=({0:f},{1:f})'.format(x[0,0],x[1,0]))

# contour plot of the cost function
plt.figure(figsize=(5,6))
delta = 0.025
x = np.arange(-1.0, 1.2, delta)
y = np.arange(-1.0, 2.0, delta)
X, Y = np.meshgrid(x, y)
c=np.array([-0.65,  -0.5, -0.3,  -0.1, 0.1,  0.3, 0.5, 0.7, 0.9, 1.1])

N=x.size
M=y.size
Z=np.zeros((M,N))

for i in range(0,M):
    for j in range(0,N):

        Z[i,j] = 0.5*(X[i,j]**2*A[0,0]+(A[0,1]+A[1,0])*X[i,j]*Y[i,j]
                 +A[1,1]*Y[i,j]**2) - X[i,j]*b[0]-Y[i,j]*b[1]

CS = plt.contour(X, Y, Z, c)
plt.clabel(CS, inline=1, fontsize=10)
plt.xlim(-1.0,1.2)
plt.ylim(-1.0,2.0)
plt.axes().set_aspect('equal', 'datalim')

# plot the trajectory
plt.plot(u[0:n+1],v[0:n+1],'-r',linewidth=2)
plt.xlabel(r'$x_1$',fontsize=14)
plt.ylabel(r'$x_2$',fontsize=14)
plt.title('Conjugate Gradient Minimization')

plt.show()