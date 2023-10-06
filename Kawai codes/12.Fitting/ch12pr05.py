#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 12.5                                                       *
%*     filename: ch12pr05.py                                              *
%*     program listing number: 12.5                                       *
%*                                                                        *
%*     This program interpolates 11-point data with the linear            *
%*     regression method.                                                 *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# Data set (no error bar)

x=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.]
y=[0.1,0.90,1.7,3.4,4.5,4.7,6.2,7.6,7.85,9.03,9.6]
x=np.array(x)
y=np.array(y)
N=y.size

# Linear regression
F=y.sum()
X=x.sum()
X2=(x**2).sum()
XF=(x*y).sum()
b=(F*X2-X*XF)/(N*X2-X**2)
a=(N*XF-X*F)/(N*X2-X**2)

# fitted curve
f=a*x+b

plt.figure(figsize=(6,5))
plt.plot(x,f,'-r')
plt.plot(x,y,'ok');
plt.xlabel('x',fontsize=14)
plt.ylabel('f(x)',fontsize=14)
plt.show()
