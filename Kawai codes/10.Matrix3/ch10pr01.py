#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 10.1 -10.3                                                 *
%*     filename: ch10pr01.m                                               *
%*     program listing number: 10.1                                       *
%*                                                                        *
%*     This program finds eigenvalues and eigenvectos of a 3x3symmetric   *
%*     matrix using the power method.                                     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/08/2013.                                    *
%**************************************************************************
"""
import numpy as np

def evpower(A):

    # tolerance
    tol=1.0e-7
    found=False

    # create arrays
    u0=np.matrix(np.zeros(3)).transpose()
    u1=np.matrix(np.zeros(3)).transpose()

    # initial guess
    x=np.matrix(np.random.rand(3)).transpose()

    # normalization
    u0=x/np.sqrt(np.asscalar(x.transpose()*x))

    # power method iteration
    n=0
    while not(found):
        x=A*x   # update x
        u1=x/np.sqrt(np.asscalar(x.transpose()*x))   # normalization
        err=np.sqrt( np.asscalar((u1-u0).transpose()*(u1-u0)) )   #error
        if err < tol:
            found=True
        else:
            u0[:]=u1[:]

        n+=1
        
    eig=np.asscalar(u1.transpose()*A*u1) 
    
    return [n, eig, u1]
        
A=np.matrix([[2.,1.,0.],[1.,2.,1.],[0.,1.,2.]])

# SOlution by eignvalue solver in numpy
eig_np, u_np = np.linalg.eig(A)

# Example 10.1
[n, eig, u] = evpower(A)
# Eigenvalue

print('\nLargest Eigenvalue (Example 10.1)')
print('Iteration=',n)
print('Eigenvale={0:f} (numpy: {1:f})'.format(eig, eig_np[0]))
print('Eigenvector=[{0:8.4f}, {1:8.4f}, {2:8.4f}]'\
                    .format(u[0,0],u[1,0],u[2,0]))
print('    (numpy):[{0:8.4f}, {1:8.4f}, {2:8.4f}]'\
      .format(u_np[0,0],u_np[1,0],u_np[2,0]))
eig_max=eig

# Example 10.2
A=np.matrix([[2.,1.,0.],[1.,2.,1.],[0.,1.,2.]])
A=np.linalg.inv(A)
[n, eig, u] = evpower(A)
# Eigenvalue
eig=1.0/eig

print('\nSmallest Eigenvalue (Example 10.2)')
print('Iteration=',n)
print('Eigenvale={0:f} (numpy: {1:f})'.format(eig, eig_np[2]))
print('Eigenvector=[{0:8.4f}, {1:8.4f}, {2:8.4f}]'\
                    .format(u[0,0],u[1,0],u[2,0]))
print('    (numpy):[{0:8.4f}, {1:8.4f}, {2:8.4f}]'\
      .format(u_np[0,2],u_np[1,2],u_np[2,2]))
eig_min=eig

# Example 10.3
A=np.matrix([[2.,1.,0.],[1.,2.,1.],[0.,1.,2.]])
eig0=0.5*(eig_min+eig_max)

A=np.linalg.inv(A-eig0*np.identity(3))
[n, eig, u] = evpower(A)
eig=1.0/eig+eig0

print('\nThe Eigenvalue between the smallest and largest (Example 10.3)')
print('Iteration=',n)
print('Eigenvale={0:f} (numpy: {1:f})'.format(eig, eig_np[1]))
print('Eigenvector=[{0:8.4f}, {1:8.4f}, {2:8.4f}]'\
                    .format(u[0,0],u[1,0],u[2,0]))
print('    (numpy):[{0:8.4f}, {1:8.4f}, {2:8.4f}]'\
      .format(u_np[0,1],u_np[1,1],u_np[2,1]))
