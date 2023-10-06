#!/usr/bin/env python
# coding: utf-8

# # Day 6 - The function

# In[1]:


def myFunction():
    print('this is a function')


# In[2]:


myFunction()


# In[24]:


def factorial(n):
    if n<0:
        print('undefined')
    else:
        f = 1
        for i in range(f,n+1):
            f = f*i
        return(f)


# In[25]:


factorial(-4)


# In[31]:


def distance(x,y,z):return (x**2+y**2+z**2)**(1/2)


# In[32]:


distance(1,1,1)


# In[42]:


def sphxyz(r, theta, phi):
    from math import sin,cos
    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)
    return (x, y, z)


# In[53]:


import math

sphxyz(1.7320508075688772,math.pi/3,math.pi/4)


# In[40]:


math.pi


# In[50]:


def xyzsph(x,y,z):
    r = (x**2+y**2+z**2)**(.5)
    theta = math.atan((x**2+y**2)**(.5)/z)
    phi = math.atan(y/x)
    return (r,theta,phi)


# In[51]:


xyzsph(1,1,1)


# In[54]:


math.pi/


# In[ ]:




