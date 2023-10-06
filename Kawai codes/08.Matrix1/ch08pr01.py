#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  8.1                                                       *
%*     filename: ch08pr01.m                                               *
%*     program listing number: 8.1                                        *
%*                                                                        *
%*     This program checks the properties of triangular matrices.         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/31/2015.                                    *
%**************************************************************************
"""
import numpy as np

# define matrices A and B (do not use array)
A=np.matrix([[2, 0, 0],[-1,1,0],[3,2,-1]])
B=np.matrix([[1, 0, 0],[2,4,0],[-1,-2,3]])

print("A*B")
print(A*B)

print('\nInverse of A')
print(np.linalg.inv(A))

print("\nDeterminant of A={0:7.5f}".format(np.linalg.det(A)))


