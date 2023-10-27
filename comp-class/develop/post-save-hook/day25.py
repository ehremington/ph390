#!/usr/bin/env python
# coding: utf-8

# # Day 25

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


# Runge-Kutta 4nd order

def f(r,t):
    x = r[0]
    y = r[1]
    fx = y
    fy = -9.81/1*np.sin(x)
    return(np.array([fx,fy],float))

# define boundary conditions

a = 0.0 # starting point
b = 10.0 # ending point
N = 1000 # number of points between a and b
dt = (b-a)/N

r = np.array([1,1],float) # initial condition

tpoints = np.arange(a, b, dt)
xpoints = []
ypoints = []

for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    k1 = dt*f(r,t)
    k2 = dt*f(r+0.5*k1,t+0.5*dt)
    k3 = dt*f(r+0.5*k2,t+0.5*dt)
    k4 = dt*f(r+k3, t+dt)
    r = r + (k1+2*k2+2*k3+k4)/6


# In[3]:


fig0, ax0 = plt.subplots(figsize=(3,3))

ax0.plot(tpoints, xpoints)
# ax0.plot(tpoints, ypoints)
ax0.set_xlabel('time')
ax0.set_ylabel(r'$\theta$')


# In[ ]:





# In[ ]:




