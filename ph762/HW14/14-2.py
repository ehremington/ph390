#!/usr/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 11:11:49 2014

@author: eric
"""
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

n = 100000
t1,t2 = 5*rand(2,n) #two lists from 0 to 5
t = np.concatenate([-t1,t2]) #this combines the two lists into one from -5 to 5
abst = abs(t)
x = np.exp(-abst**2)   #forms a normal distribution from this list of values

mean = np.sum(x)/(2*n)
print(mean)
