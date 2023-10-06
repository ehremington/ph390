# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example  9.1                                                       *
%*     filename: ch09pr01.py                                              *
%*     program listing number: 9.1                                        *
%*                                                                        *
%*     This program solves a two-dimensional non-linear equation          *
%*     by newton-raphson method.                                          *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/10/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# set a tolerance
tol = 1.0e-8

# step factor (between 0 and 1)
a = 0.1

#define functions
def f(x,y):
    return [2.0*x + 3.0*x*y - 1.0, x*y + 3.0*y + 1.0]

#define Jacobian
def J(x,y):
    return [[2.0+3.0*y,3.0*x],[y,x+3.0]]

nmax=1000
x=np.zeros(nmax+1)
y=np.zeros(nmax+1)
b=np.zeros(2)
A=np.zeros((2,2))

# iniital guess
x[0] = 1.0
y[0] = 0.0

b=-np.array(f(x[0],y[0]))

err=np.sqrt(b[0]*b[0]+b[1]*b[1])
if err < tol:
    found = True
else:
    found = False

n=0
while not(found) and n<nmax:
    
    # Construct linear equation
    A = np.array(J(x[n],y[n]))

    # solve the linear equation
    z=np.linalg.solve(A,b)
    x[n+1]=x[n]+a*z[0]
    y[n+1]=y[n]+a*z[1]
    
    # check error
    b=-np.array(f(x[n+1],y[n+1]))
    err=np.dot(b,b)
    print('err=',err)
    if err < tol:
       found = True
    else:
       n+=1


plt.figure(figsize=(6,5))     
plt.plot(x[0:n],y[0:n],'ok',linewidth=2)
plt.axes().set_aspect('equal', 'datalim')
plt.text(x[0]+0.01,y[0],'Start Here')
plt.text(x[n]+0.01,y[n]-0.01,'Converged',color='r')
plt.xlabel('x',fontsize=14)
plt.ylabel('y',fontsize=14)
plt.show()

xa = (-1.0+np.sqrt(7.0))/2.0
ya = (-5.0+np.sqrt(7.0))/9.0
print('Iternations = {0:d}'.format(n))
print('Solution= {0:f}, {1:f}'.format(x[n],y[n]))
print('   Exact= {0:f}, {1:f}'.format(xa,ya))
