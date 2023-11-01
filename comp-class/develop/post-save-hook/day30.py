#!/usr/bin/env python
# coding: utf-8

# # Day 30

# In[4]:


import numpy as np
import matplotlib.pyplot as plt


# In[6]:


fig0,ax0 = plt.subplots()
x = np.linspace(0,1,100)
ax0.plot(x, np.exp(-x))
ax0.plot(x,x)


# In[8]:


x = 1
for i in range(20):
    x = np.exp(-x)
    print(x)


# In[11]:


x = 100
for i in range(20):
    x = np.exp(-x)
    print(x)


# Let's try with a new function:
# $$x = e^{1-x^2}$$

# In[15]:


fig2, ax2 = plt.subplots()
x = np.linspace(0,2,100)
ax2.plot(x,np.exp(1-x**2))
ax2.plot(x,x)


# In[14]:


x = 0.5
for i in range(20):
    x = np.exp(1-x**2)
    print(x)


# Let's try the inverted form of the previous equation:
# $$x=\sqrt{1-\ln x}$$

# In[17]:


fig1, ax1 = plt.subplots()
x = np.linspace(0,2,100)
ax1.plot(x, np.sqrt(1-np.log(x)))
ax1.plot(x,x)


# In[18]:


x = 0.5
for i in range(20):
    x = np.sqrt(1-np.log(x))
    print(x)


# In[21]:


x = 2 # guess
diff = 1 # arbitrary
while abs(diff) > 2**-32:
    x1 = np.sqrt(1-np.log(x))
    diff = x - x1
    x = x1
    
print(x)


# In[2]:


def bisection(function, lower_guess, upper_guess, tolerance=2**-32):
    midpoint = (lower_guess + upper_guess)/2
    while upper_guess - lower_guess > tolerance:
        if function(lower_guess)*function(midpoint)<0:
            upper_guess = midpoint
            midpoint = (lower_guess + upper_guess)/2
        elif function(midpoint)*function(upper_guess)<0:
            lower_guess = midpoint
            midpoint = (lower_guess + upper_guess)/2
        elif function(lower_guess)*function(midpoint)>0 and function(midpoint)*function(upper_guess)>0:
            print('no unique root in that bracket')
            break
    return(midpoint)


# In[3]:


def func(x):
    return(x-np.exp(1-x**2))

bisection(func, 0, 1.7)


# In[33]:


func(2)


# In[ ]:




