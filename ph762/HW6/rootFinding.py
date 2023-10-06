#!/usr/bin/python3

#This is a module for root finding and containts two methods

def x_mid(lower_x,upper_x): return (lower_x+upper_x)/2.0

def bisection(function, lower_x, upper_x,tolerance = 2**-32):
    mid_x = x_mid(lower_x,upper_x)
    bracket = [lower_x,mid_x,upper_x]
    y_bracket = [function(bracket[0]),function(bracket[1]),function(bracket[2])]
    n = 0
    r = bracket[1]
    while abs(function(bracket[1]))>tolerance:
        if y_bracket[0]*y_bracket[1]<0:
            x_out = bracket.pop(2)
            x_in = x_mid(bracket[0],bracket[1])
            bracket.insert(1,x_in)    
        elif y_bracket[1]*y_bracket[2]<0:
            x_out = bracket.pop(0)
            x_in = x_mid(bracket[0],bracket[1])
            bracket.insert(1,x_in)
        elif y_bracket[0]*y_bracket[1]>0 and y_bracket[1]*y_bracket[2]>0:
            print("no unique root in that bracket")
            break
        y_bracket = [function(bracket[0]),function(bracket[1]),function(bracket[2])]
        #print(bracket,y_bracket,n)
        n = n+1
        r = bracket[1]   #root is r
    return r, n, [bracket[1],bracket[2]] 






def newtonRaphson(f, df, guess, tolerance = 2**-32):
    x = guess
    n = 0
    while abs(f(x))>tolerance:
        x = x - f(x)/df(x)
        n = n+1
    root = x
    return root,n

def secant(function,guess,delta, tolerance = 2**-32):
    x0 = guess
    x1 = x0 + delta
    n = 0
    while abs(function(x1))>tolerance:
        xt = x1
        #print(x0,x1,function(x0),function(x1))
        x1 = x1 - (x1-x0)/(function(x1)-function(x0))*function(x1)
        x0 = xt
        n = n+1
    root = x1
    return root,n