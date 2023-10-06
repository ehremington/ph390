#!/usr/bin/python3
"""
Created on Thu Mar  6 15:52:57 2014

This contains all derivative taking functions of one variable

@author: eric
"""


def back(function,x_value, stepSize):
    return (function(x_value)-function(x_value-stepSize))/stepSize

def forward(function,x_value, stepSize):
    return (function(x_value+stepSize)-function(x_value))/stepSize

def mid(function,x_value, stepSize):
    return (function(x_value+stepSize)-function(x_value-stepSize))/(2*stepSize)

def doubled(function, x_value, stepSize):
    return (function(x_value+stepSize)-2*function(x_value)+function(x_value-stepSize))/(stepSize**2)