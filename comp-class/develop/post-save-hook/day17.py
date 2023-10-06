#!/usr/bin/env python
# coding: utf-8

# # Day 17

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


def integrate(func, a, b, steps):
    # using trapezoid method
    h = (b-a)/steps
    s = (func(a)+func(b))*h/2
    x = np.linspace(a, b, steps+1)
    y = func(x)
    s = s + np.sum(y[1:steps]*h)
    return(s)

def integrate1(func, a, b, steps):
    if steps%2 != 0:
        steps = steps + 1
    h = (b-a)/steps
    x = np.linspace(a, b, steps+1)
    y = func(x)
    s = h/3*(y[0] + y[-1] + 4*np.sum(y[1:steps:2]) +
             2*np.sum(y[2:steps-1:2]))
    return(s)


# In[3]:


def f(x):
    return(x**4-2*x+1)


# In[4]:


integrate1(f, 0, 2, 101)


# In[5]:


def f1(t):
    return(np.exp(-(t**2)))

def capitalE(x):
    return(integrate1(f1,0,x,100))


# In[6]:


fig0, ax0 = plt.subplots()

x = np.linspace(0, 3, 31)
y = capitalE(x)

ax0.plot(x,f1(x), 'o', label=r'f1(t)')
#ax0.plot(x,y, 'o', label=r'$\int f1(t) dt$')


# In[7]:


capitalE(3.3)


# In[8]:


y[30]


# In[9]:


x


# In[10]:


capitalE(x)


# In[ ]:




