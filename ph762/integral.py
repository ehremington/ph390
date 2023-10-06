#!/usr/bin/python3
"""
Created on Thu Mar  6 16:19:53 2014

This contains all integrating functions of one variable

@author: eric
"""
import numpy as np
   
def rect(function,start,stop,steps):
    h = (stop-start)/steps
    x = np.linspace(start,stop,steps)
    f = function(x)
    rect = np.sum(f)*h
    return rect
    
def trap(function,start,stop,steps):
    h = (stop-start)/steps
    x = np.linspace(start,stop,steps)
    f = function(x)
    rect = np.sum(f[1:steps-1])*h
    trap = rect + h*(f[0]+f[steps-1])/2
    return trap
    
def simpson(function,start,stop,steps):
    if steps%2 == 0:
        steps = steps + 1   
    h = (stop-start)/steps    
    x = np.linspace(start,stop,steps)
    f = function(x)
    simp = (f[0]+f[steps-1]+5*(f[1]+f[steps-2])+6*np.sum(f[2:steps-2]))*h/3
    return simp

def GaussLagQuad8(function):
    x=np.asarray([1.7027963230510100e-1, 9.0370177679937991e-1,\
       2.2510866298661307,    4.2667001702876588,\
       7.0459054023934657,    1.0758516010180995e+1,\
       1.5740678641278005e+1, 2.2863131736889264e+1])
    w=np.asarray([3.6918858934163753e-1, 4.1878678081434296e-1,\
       1.7579498663717181e-1, 3.3343492261215652e-2,\
       2.7945362352256725e-3, 9.0765087733582131e-5,\
       8.4857467162725315e-7, 1.0480011748715104e-9])
    integral = np.sum(w*function(x))
    return(integral)

def gaussHermQuad8(function):
    x = np.asarray([-0.38118699,-1.157193712,-1.981656757,-2.93063742,0.38118699,1.157193712,1.981656757,2.93063742])
    w = np.asarray([0.661147013,0.207802326,0.017077983,0.000199604,0.661147013,0.207802326,0.017077983,.000199604])
    integral = np.sum(w*function(x))
    return(integral)

def gaussLegQuad8(function):
    x = np.asarray([-0.183434643,-0.52553241,-0.796666477,-0.960289857,0.183434643,0.52553241,0.796666477,0.960289857])
    w = np.asarray([0.362683783,0.313706646,0.222381035,0.101228536,0.362683783,0.313706646,0.222381035,0.101228536])
    integral = np.sum(w*function(x))




'''
###########____depricated code below_____####################

def listFromFuncEndpoints(function,x_begin,x_end,steps):    
    x_list = np.linspace(x_begin,x_end,steps)
    y_list = function(x_list)
    return x_list, y_list


def TrapFromList(x_list,y_list):
    integral = 0.0
    for i in range(len(x_list)-1):
        integral = integral + (y_list[i+1]+y_list[i])*(x_list[i+1]-x_list[i])/2.0
    return integral

def TrapFromFunction(function,start,end,steps):
    x_list, y_list = listFromFuncEndpoints(function,start,end,steps)
    integral = IntTrapFromList(x_list,y_list)
    return integral

def SimpFromList(x_list, y_list):
    if len(x_list)%2 == 0:
        print('list does not have an odd number of elements')
    else:
        integral = 0.0
        for i in range(1,len(x_list)-1,2):
            integral = integral + (y_list[i-1]+4*y_list[i]+y_list[i+1])*(x_list[i+1]-x_list[i-1])/6.0        
        return integral

def SimpFromFunction(function,start,end,steps):
    #checks for an even number of steps
    if steps%2 == 0:
        steps = steps + 1
        #print('I fixed incorrect number of steps by adding one more')
    x_list, y_list = listFromFuncEndpoints(function,start,end,steps)
    integral = IntSimpFromList(x_list,y_list)
    return integral
    
'''