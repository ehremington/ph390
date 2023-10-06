#!/usr/bin/python3

import rootFinding as root
import math
import matplotlib.pyplot as plt
import numpy as np

#define a function to find the root
def f(z):
    return z/np.tan(z)+(6**2-z**2)**(2**-1)
def df(z):
    return 1/np.tan(z) + z/(np.sin(z)**2) - (z**2)*(6**2-z**2)**(-2**-1) 



#this section plots the function so we can properly bracket the problem

x = np.arange(0,6,.1)
y = f(x)
zero = np.zeros(np.shape(x))
plt.plot(x,y,x,zero)
plt.show()



   
#this does about 12 steps of bisection and then passes off the final bracket to secant method    
print('The first root is bracketed bases on inspection of the graph')
r1,n1,bracket1 = root.bisection(f,2.5,3,2**-6)
print('12ish steps of bisection method','\n',r1,'(iterations = %i)'%n1)
#this is a bad name but this is just a refinement of the first root
r11,n11 = root.secant(f,bracket1[0],(bracket1[1]-bracket1[0])/10,2**-32)
print('further refinement by secant method','\n',r11,'(iterations = %i)'%n11)
#this is how long the bisection method would have taken by itself
r111,n111,bracket111 = root.bisection(f,2.5,3,2**-32)
print('how long bisection method would take by itself for this precision:',n111)

print('\n Now for the second root, bracketed by inspection of the graph')
#the next root, estimated by the graph
r2,n2,bracket2 = root.bisection(f,5,5.5,2**-6)
print('12ish steps of bisection method','\n',r2,'(iterations = %i)'%n2)
#this is a bad name but this is just a refinement of the first root
r22,n22 = root.secant(f,bracket2[0],(bracket2[1]-bracket2[0])/10,2**-32)
print('further refinement by secant method','\n',r22,'(iterations = %i)'%n22)
#this is how long the bisection method would have taken by itself
r222,n222,bracket222 = root.bisection(f,5,5.5,2**-32)
print('how long bisection method would take by itself for this precision:',n222)
