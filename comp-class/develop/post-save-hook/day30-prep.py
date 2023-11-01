#!/usr/bin/env python
# coding: utf-8

# ## Day 30 prep

# In[17]:


get_ipython().run_line_magic('matplotlib', 'widget')
import numpy as np
import matplotlib.pyplot as plt


# In[30]:


fig1, ax1 = plt.subplots()
x = np.linspace(0,2,100)
ax1.plot(x,np.exp(-x))
ax1.plot(x,x)


# In[18]:


x = 1
for i in range(10):
    x = np.exp(-x)
    print(x)


# In[19]:


x = -20
for i in range(20):
    x = np.exp(-x)
    print(x)


# Let's try a slightly different function though.
# $$x = e^{1-x^2}$$

# In[28]:


x = 10
for i in range(20):
    x = np.exp(1-x**2)
    print(x)


# The above function doesn't seem so different, and yet the results aren't working. Let's plot and see.

# In[34]:


fig0, ax0 = plt.subplots()
t = np.linspace(0,2,100)
ax0.plot(t, np.exp(1-t**2))
ax0.plot(t,t)


# Now the problem relaxation method doesn't seem to work so well. However, there is an interesting trick with this we can play. Let's solve the problem the other way around.
# 
# $$x = e^{1-x^2} \rightarrow x = \sqrt{1-\log{x}}$$

# In[55]:


fig2, ax2 = plt.subplots()
x = np.linspace(0.01,2,100)
ax2.plot(x, np.sqrt(1-np.log(x)))
ax2.plot(x,x)


# In[38]:


x = 2
for i in range(20):
    x = np.sqrt(1-np.log(x))
    print(x)


# We can beef up this method slightly to report a result after the answer get sufficiently close to the "true" answer. We can do this by subtraction for example.

# In[50]:


x = 2
e = 1
while abs(e) > 0.0000001:
    x1 = np.sqrt(1-np.log(x))
    e = x1-x
    x = x1
print(x)


# In[47]:


x+e


# In[ ]:




