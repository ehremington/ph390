#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 10.4                                                       *
%*     filename: ch10pr02.py                                              *
%*     program listing number: 10.2                                       *
%*                                                                        *
%*     This program finds eigenvalues and eigenvectos of a symmetric      *
%*     matrix using the inverse iteration method method.                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/17/2017.                                    *
%**************************************************************************
"""
import numpy as np

# Set a matrix
A=np.matrix([[2,1,0],[1,2,1],[0,1,2]])

for q in (0.5,1.0,1.5,2.5,3.0,3.5):
    
    #tolerance
    tol = 1e-7

    found=False
    n=1

    # Generate a random vector
    b0=np.random.rand(3,1)
    b0=b0/np.linalg.norm(b0)

    B=A-q*np.identity(3)

    # inverse iteration method
    while not(found):
        y=np.linalg.solve(B,b0)
        b1=y/np.linalg.norm(y)
        if b1[0]<0: # correct the phase.
            b1=-b1

        err=np.linalg.norm(b1-b0)
        if err < tol:
            found=True
        else:
            b0=b1
            n+=1

    # eigenvalue
    eig= np.asscalar(b1.transpose()*A*b1)
    print('Guess={0:6.3f},   Eigenvalue={1:10.6f}'.format(q,eig))

