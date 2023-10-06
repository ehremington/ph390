#!/usr/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 20:00:57 2014

@author: eric
"""
import numpy as np
from numpy.random import rand


n = 100000
x = rand(n)
mean = np.sum(x)/n
variance = np.sum((x-mean)**2)/n
print(mean,variance)
print('theoretical mean: 1/2, theoretical variance: 1/12')