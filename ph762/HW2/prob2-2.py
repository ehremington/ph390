#!/usr/bin/python
#import modules
import calculus as calc
import matplotlib.pyplot as plot
import numpy as np

print "This program evaluates the second derivative of 1/12*x**4 at any position"
print "It also prints a table of these values at various step sizes"
print "Not only that, but it will graph the error for comparison"
#need the function to evaluate
def func(x):
    return 1/12.0*x**4


#get position
x=1

h = []
diffTable = []
error = []

print "h               " + "derivative           " + "error"
for i in range(21):
    h.append(10**(-i))
    diffTable.append(calc.doubleDiff(func,x,h[i]))
    error.append(abs(1.0-diffTable[i]))
    print "%.2e         "%h[i] + "%f           "%diffTable[i] + "%f" %error[i] 
    
plot.plot(h,error, 'o')
plot.xscale('log')
plot.yscale('log')
plot.axis([10**-21,10,10**-16,100])
plot.show()

