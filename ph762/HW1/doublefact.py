#!/usr/bin/python
#########################################################################
#
#		A program to calculate n!!
#		Eric Remington	
#		Computational Physics PH 762
#		1/14/14	
#
#########################################################################

import math

#user input
n = eval(raw_input("What integer should I double factorial? "))


#determine size of n to decide approximation sceme

if n<= 30:
    #we work with the log of the input value, so just the exponents
    exponent = math.log10(n)
	#this sums the logs of the factorial terms to obtain an approximation
    while n>2:
        exponent = exponent+math.log10(n-2)
        n = n-2

    expToDispl = math.floor(exponent)
    mantissa = 10**(exponent-expToDispl)
    print "%.3f x 10^%i" % (mantissa, expToDispl)
    print "done"

#even values computed here
elif n%2==0:
    exponent = (n*.5)*math.log(2.0) + .5*math.log(n*math.pi) + n*.5*math.log(n*.5) - n*.5 + math.log(1 + 1/(6*n))
    #change of base needed here
    base10exp = exponent/math.log(10)
    #dislay mantissa and exponent in a nice form
    expToDispl = math.floor(base10exp)
    mantissa = 10**(base10exp-expToDispl)
    print "%.3f x 10^%i" % (mantissa, expToDispl)
    print "done"
 
#odd values computed here but otherwise same as even case    
else:
    exponent = (1-n)*.5*math.log(n-1)+n*math.log(n)-(n+1)*.5
    base10exp = exponent/math.log(10)

    expToDispl = math.floor(base10exp)
    mantissa = 10**(base10exp-expToDispl)
    print "%.3f x 10^%i" % (mantissa, expToDispl)
    print "done"
    
