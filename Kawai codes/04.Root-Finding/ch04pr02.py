#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*  Example  4.2                                                          *
%*  filename: ch04pr02.py                                                 *
%*  program listing number: 4.2                                           *
%*                                                                        *
%*     This program finds roots of a cubic equation                       *
%*               a*x^3 + b*x^2 + c*x + d = 0                              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/14/2017.                                    *
%**************************************************************************
"""
import numpy as np
# coefficients
b=-9.0; c=23.0; d=-15.0
k=0

# formula for cubic polynomials
F=(3.0*c-b**2)/3.0
G=(2.0*b**3 - 9.0*b*c + 27.0*d)/27.0
H=G**2/4.0 + F**3/27.0

yes = False

if H>0.0:
    S= (np.sqrt(H)-G/2.0)**(1./3.)
    U=-(np.sqrt(H)+G/2.0)**(1./3.)
    x1 = S+U-b/3.0
    print("Answer = {0:10.5f}".format(x1))

else:
    I=np.sqrt(G**2/4.0-H)
    J=I**(1./3.)
    K=np.arccos(-G/(2.0*I))
    M=np.cos(K/3.0)
    N=np.sqrt(3)*np.sin(K/3.0)
    P=-b/3.0
    x3=np.zeros(3)    
    x3[0] = P-J*(M+N)
    x3[1] = P-J*(M-N)
    x3[2] = P+2.0*J*M
    print("Answer = {0:10.5f}, {1:10.5f}, {2:10.5f}"
          .format(x3[0],x3[1],x3[2]))

  