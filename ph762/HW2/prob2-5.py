#! /usr/bin/python
"""
This program calculates the heat capacity of a free electron gas using
normalized coordinated for heat capacity and temperature.

The normalized heat capacity will be plotted against normalized temperature
from Tn = 0 to Tn = 1

The numpy python package is required.
"""
import numpy as np
import calculus as cal
import matplotlib.pyplot as plt
import math
"""This integration will take place in two parts. Values for Tn near 1 will be
evaluated with a combination of Gaussian-Laguere quadrature and Simpson's Rule
and for Tn near 0, by a substitution trick and Simpson's rule.

A list from the integration routines will be populated and graphed at the end
"""

#this makes a numpy array
Tn = np.linspace(.01,1.0,1000)
#Notice the above starts at .01 instead of 0. This is due to the smallness of these functions at values less than that


#This function is the integrand for the small values of Tn
def func1(s):
    return (1/(s**4))*math.e**(-1/s)/(math.e**(-1/s)+1)**2

#This function call the correct integration method and inputs the Tn as limits of integration
def smallTnIntegral(tn):
    return cal.IntSimpFromFunction(func1, .01, tn, 10001)
"""
#This performs the integration and multiplication and adds the results to another array
Ce = Tn*math.pi**2/3 - Tn*smallTnIntegral(Tn)
"""
Ce =[]
#need to redo this with a for loop to get this to work
for tn in Tn:
    Ce.append(tn*math.pi**2/3 - tn*smallTnIntegral(tn))
    


#plotting performed here
plt.plot(Tn,Ce)
plt.xlabel("normalized Temp")
plt.ylabel('normalized C')
plt.show()

if __name__ == "__main__":
    print("done")