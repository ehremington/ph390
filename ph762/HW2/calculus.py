import numpy as np


def backDiff(function,x_value, stepSize):
    return (function(x_value)-function(x_value-stepSize))/stepSize

def forwardDiff(function,x_value, stepSize):
    return (function(x_value+stepSize)-function(x_value))/stepSize

def midDiff(function,x_value, stepSize):
    return (function(x_value+stepSize)-function(x_value-stepSize))/(2*stepSize)

#returns odd number lists of a function some steps size around a central value
def listFromFuncCentralpoint(function,centralValue, steps, stepSize):
    x_list = []
    y_list = []
    for i in range(2*steps+1):
        x_list.append(centralValue+(i-steps)*stepSize)
        y_list.append(function(centralValue+(i-steps)*stepSize))
        
    return x_list, y_list

def listFromFuncEndpoints(function,x_begin,x_end,steps):    
    x_list = np.linspace(x_begin,x_end,steps)
    y_list = function(x_list)
    # old code i don't want to get rid of below
    """
    stepSize = (x_end-x_begin)/float(steps)
    for i in range(steps+1):
        x_list.append(x_begin+i*stepSize)
        y_list.append(function(x_begin+i*stepSize))"""
    return x_list, y_list
    
        
def doubleDiff(function, x_value, stepSize):
    x_list,y_list = listFromFuncCentralpoint(function,x_value,2,stepSize)
    y_prime = []
    for i in range(1,4):
        y_prime.append((y_list[i+1]-y_list[i-1])/(2*stepSize))
    y_doubleprime = (y_prime[2]-y_prime[0])/(2*stepSize)
    return y_doubleprime

def IntTrapFromList(x_list,y_list):
    integral = 0.0
    for i in range(len(x_list)-1):
        integral = integral + (y_list[i+1]+y_list[i])*(x_list[i+1]-x_list[i])/2.0
    return integral

def IntTrapFromFunction(function,start,end,steps):
    x_list, y_list = listFromFuncEndpoints(function,start,end,steps)
    integral = IntTrapFromList(x_list,y_list)
    return integral

def IntSimpFromList(x_list, y_list):
    if len(x_list)%2 == 0:
        print('list does not have an odd number of elements')
    else:
        integral = 0.0
        for i in range(1,len(x_list)-1,2):
            integral = integral + (y_list[i-1]+4*y_list[i]+y_list[i+1])*(x_list[i+1]-x_list[i-1])/6.0        
        return integral

def IntSimpFromFunction(function,start,end,steps):
    #checks for an even number of steps
    if steps%2 == 0:
        steps = steps + 1
        print('I fixed incorrect number of steps by adding one more')
    
    
    x_list, y_list = listFromFuncEndpoints(function,start,end,steps)
    integral = IntSimpFromList(x_list,y_list)
    return integral
    
def IntGaussLagQuad8(function):
    x=np.asarray([1.7027963230510100e-1, 9.0370177679937991e-1,\
       2.2510866298661307,    4.2667001702876588,\
       7.0459054023934657,    1.0758516010180995e+1,\
       1.5740678641278005e+1, 2.2863131736889264e+1])
    w=np.asarray([3.6918858934163753e-1, 4.1878678081434296e-1,\
       1.7579498663717181e-1, 3.3343492261215652e-2,\
       2.7945362352256725e-3, 9.0765087733582131e-5,\
       8.4857467162725315e-7, 1.0480011748715104e-9])
    integral = sum(w*function(x))
    return integral