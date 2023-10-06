#!/usr/bin/python

#########################################################################
#
#		A program to compare different methods for solving the
#			quadratic equation
#		Eric Remington	
#		Computational Physics PH 762
#		1/19/14	
#
#########################################################################

import sys

#obtain user input for values of b and c
print("This program prints a table of solutions to the quadratic equation a*x^2 + b*x + c = 0\nfor values of a ranging from 0.1 to machine epsilon")
b = eval(input("Provide a positive value for b:"))
c = eval(input("Provide a positive value for c:"))

#initialize list
ep = []
i=0.1

#populate list with values of a from 0.1 to machine epsilon
while i > sys.float_info.epsilon:
     ep.append(i)
     i = i/10.0

#write table of solutions for the two methods
for j in ep:
	quad1 = (-b + (b**2 - 4*j*c)**.5)/(2*j)
	quad2 = -2*c/(b+(b**2-4*j*c)**.5)
	print("%e          "%j + "%f          "%quad1 +"%f"%quad2)

