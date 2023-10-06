
"""
%*********************************************************************
%*     Example  1.7                                                  *
%*     filename: ch01pr01.py                                         *
%*     program listing number: 1.1                                   *
%*                                                                   *
%*     This program finds a machine epsilon by evaluating            *
%*                                                                   *
%*           1 + 2^(-n) > 1                                          *
%*                                                                   *
%*     At a certain positive n, this inenqualty becomes false.       *
%*     Then, the machine epsilon is 2^(n-1).                         *
%*                                                                   *
%*     Programed by Ryoichi Kawai for Computational Physics Course   *
%*     Revised on 12/27/2016                                         *
%*********************************************************************
"""
# Load NumPy package
import numpy as np

# Find the machine epsilon for 64 bit float
epsilon = 1.0  # create a float64 variable
n = 0          # reset a counter

# Reduce the value of epsilon until it becomes too small
while 1.0+epsilon > 1.0:
   epsilon = epsilon/2.0
   n = n+1

# The smallest single floating value which can be added to one.
epsilon = epsilon+epsilon

# Show the results
print("Machine epsilon for 64 bit floating point")
print("Stopped after {0:3d} itersations".format(n))
print("machine epsilon by computation = {0:16.7e}".format(epsilon))
print("machine epsilon by Numpy       = {0:16.7e}".format(np.finfo(np.float).eps))
print("1+epsilon   = {0:24.20e}".format(1+epsilon))
print("1+epsilon/2 = {0:24.20e}".format(1+epsilon/2.0))
