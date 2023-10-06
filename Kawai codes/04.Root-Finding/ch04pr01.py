#!/usr/bin/env python3
"""
%**************************************************************************
%*     Example  4.1                                                       *
%*     filename: ch04pr01.py                                              *
%*     program listing number: 4.1                                        *
%*                                                                        *
%*     This program finds roots of a quadratic equation                   *
%*               a x^2 + x + 1/4 = 0                                      *
%*     for which the quadratic equation formula                           * 
%*            x = (-1 + sqrt(1 - a))/(2a)                                 *
%*     fails for small value of a.                                        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/07/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

b=1.0
c=0.25
a=1.0
x=np.zeros(50)
y1=np.zeros(50)
y2=np.zeros(50)
n=0
while a > np.finfo(float).eps:
    x[n]=a
    d=b**2-4*a*c
    y1[n]=(-b+np.sqrt(d))/(2.0*a)
    y2[n]=-2*c/(b+np.sqrt(d))
    print("a={0:22.16e}, regular={1:22.16e}, smart={2:22.16e}"
          .format(x[n],y1[n],y2[n]))
    n+=1
    a=a/10.

plt.ioff()
plt.figure(figsize=(6,5))
plt.loglog(x[0:n],abs(y2[0:n]+0.25), 'ob', label='smart')
plt.loglog(x[0:n],abs(y1[0:n]+0.25), 'or', label='regular')
plt.legend(loc=4)
plt.xlabel('a')
plt.ylabel('x')
plt.show()
