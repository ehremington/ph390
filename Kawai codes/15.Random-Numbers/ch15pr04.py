#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 15.4                                                       *
%*     filename: ch15pr04.m                                               *
%*     program listing number: 15.4                                       *
%*                                                                        *
%*     This program evaluates 1st through 4th moments of normal           *
%*     distribution.                                                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
"""
import numpy as np

N=100000

v=np.random.normal(0.0,1.0,N)
m=np.zeros(4)
print('order ','  moment  ')
for i in [1,2,3,4]:
    m[i-1]=sum(v**i)/N
    print('{0:3d}   {1: 10.4e}'.format(i,m[i-1]))
    
r=m[3]/m[1]**2
print('\n m4/m2**2={0:8.4f} (exact=3)'.format(r))

