#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt


# In[ ]:


# define differential equation

def f(x,t):
    return(-x**3 + np.sin(t))

# define boundary conditions

a = 0.0 # starting point
b = 10.0 # ending point
N = 10000 # number of points between a and b
dt = (b-a)/N

x = 0.0 # initial condition


tpoints = np.arange(a, b, dt)
xpoints = []

for t in tpoints:
    xpoints.append(x)
    x = x + dt*f(x,t)


# In[17]:


fig0, ax0 = plt.subplots()

ax0.plot(tpoints, xpoints)
ax0.set_xlabel('time')
ax0.set_ylabel('position')


# In[ ]:





# In[25]:


# define second order differential equation

def f(x,v,t):
    return(-9.8)

# define boundary conditions

a = 0.0 # starting point
b = 1.45 # ending point
N = 10000 # number of points between a and b
dt = (b-a)/N

x = 10.0 # initial condition of position
v = 0.0 # initial condition of velocity

tpoints = np.arange(a, b, dt)
xpoints = []
vpoints = [] 

for t in tpoints:
    xpoints.append(x)
    vpoints.append(v)
    v = v + dt*f(x,v,t)
    x = x + dt*v


# In[26]:


fig1, ax1 = plt.subplots()

ax1.plot(tpoints, xpoints)
ax1.set_xlabel('time')
ax1.set_ylabel('position')


# In[27]:


xpoints


# In[ ]:





# In[ ]:




