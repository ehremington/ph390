#!/usr/bin/env python
# coding: utf-8

# ## Day 21
# 
# Solving ODE's using python rather than excel

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[21]:


from math import sin
from numpy import arange
#from pylab import plot,xlabel,ylabel,show
from matplotlib.pyplot import subplots

def f(x,t):
    return -x**3 + sin(t)

a = 0.0           # Start of the interval
b = 10.0          # End of the interval
N = 100          # Number of steps
h = (b-a)/N       # Size of a single step
x = 0.0           # Initial condition

tpoints = arange(a,b,h)
xpoints = []
for t in tpoints:
    xpoints.append(x)
    x += h*f(x,t)

fig, ax = subplots()
ax.plot(tpoints,xpoints)
ax.set_xlabel("t")
ax.set_ylabel("x(t)")


# In[3]:


len(xpoints)


# In[8]:


tpoints[0:10]


# In[9]:


xpoints[0:10]


# In[ ]:




