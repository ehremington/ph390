#!/usr/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 20:00:09 2014

@author: eric
"""
import numpy as np
import matplotlib.pyplot as plt


numberOfParticles = 1000000
f_list = []
for part in range(numberOfParticles):
    steps = 2*np.random.random_integers(0,1,1000)-1
    final_position = np.sum(steps)
    f_list.append(final_position)

#plt.hist(f_list,50)
#plt.show()
mean = sum(f_list)/numberOfParticles
f_listarray = np.asarray(f_list)
variance = sum((f_listarray)**2)/numberOfParticles
skewness = sum((f_listarray)**3)/(numberOfParticles*variance**(3/2))
kurtosis = sum((f_listarray)**4)/(numberOfParticles*variance**2)
print('mean = ',mean,'variance = ',variance,'skewness = ',skewness,'kurtosis = ',kurtosis)
print('theoretical values for comparison: mean = 0, variance = 1000 (Number of steps),  skewness = 0, kurtosis = 3')