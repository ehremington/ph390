#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Tue Feb 11 13:13:17 2014

@author: eric
"""
import sys


#initialize list
ep = []
i=0.1

#populate list with values of a from 0.1 to machine epsilon
while i > sys.float_info.epsilon:
     ep.append(i)
     i = i/10.0
#this should add one more value beyond the machine epsilon
ep.append(i/10)

a = 1.0
b = 2.0

print('value of c             '+'method 1            '+'method 2')
for c in ep:
	quad1 = (-b + (b**2 - 4*a*c)**.5)/(2*a)
	quad2 = -2*c/(b+(b**2-4*a*c)**.5)
	print("%e          "%c + "%f          "%quad1 +"%f"%quad2)


print('The reason the traditional quadratic formula does not fail in this case \n \
       is due to the fact that there is an actual root at zero, so the numerator \n\
       happens to be going to zero due to round off error where there is actually a root')