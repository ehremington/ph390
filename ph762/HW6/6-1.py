#!/usr/bin/python3
# -*- coding: utf-8 -*-

from ode import numerov
from ode import cRK4finalpoint
import rootFinding as rf


t0 = -5
x0 = 0
y0 = 2**-10
tf = 0
xf = 0
yf = 0
def v(t):
    if abs(t) < 1:
        return 0
    else:
        return 36


def eigen(t_f, x_f, y_f, t_0, x_0, y_0, lambda_guess1,lambda_guess2, y_b_c=False, tolerance = 2**-25):
    """
    the initial condition list must be in the order of tf,xf,t0,x0,y0
    """
    def shot(lam):
        def s(t):return 0
        def w(t):return -(v(t)-lam)
        def f(x,y,t): return y
        def g(x,y,t): return s(t)-w(t)*x
        if y_b_c: #tests what kind of boundary condition to shoot for
            xf_test,yf_test, tf_test = cRK4finalpoint(f,g,t_f,y_0,x_0,t_0,2**-10)
            accuracy = y_f-yf_test
        else:
            xf_test,yf_test, tf_test = cRK4finalpoint(f,g,t_f,y_0,x_0,t_0,2**-10)
            accuracy = x_f-xf_test
        return accuracy
    lam_temp,n_temp,bracket = rf.bisection(shot,lambda_guess1,lambda_guess2,2**-3)
    lam_final,n = rf.secant(shot,lam_temp,(bracket[1]-bracket[0])/10,tolerance)
    #(bracket[1]-bracket[0])/10
    return lam_final, n    
    #return lam_temp, n_temp#, bracket
print('wait for it...')    
e1,n1 = eigen(tf,xf,yf,t0,x0,y0,1,2,True)
e2,n2 = eigen(tf,xf,yf,t0,x0,y0,7,8)
e3,n3 = eigen(tf,xf,yf,t0,x0,y0,12,20,True)
e4,n4 = eigen(tf,xf,yf,t0,x0,y0,20,30)
print(e1,e2,e3,e4)
