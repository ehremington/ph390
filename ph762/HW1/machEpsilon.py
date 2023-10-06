#!/usr/bin/python
#########################################################################
#
#		A program to calculate machine epsilon
#		Eric Remington	
#		Computational Physics PH 762
#		1/9/14	
#
#########################################################################


#define a starting point for machine epsilon
epsilon = 1.
#in order to know how many iterations it takes to get to machine epsilon from unity
n=0
while 1+epsilon/2>1:
    epsilon = epsilon/2
    n = n+1

 
print epsilon
print n
