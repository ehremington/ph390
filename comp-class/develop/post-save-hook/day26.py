#!/usr/bin/env python
# coding: utf-8

# # Day 26

# In[ ]:





# In[5]:


# Runge-Kutta 4nd order

def f(r,t):
    y = r[0]
    v = r[1]
    fy = v
    fv = -9.81
    return(np.array([fy,fv],float))

# define boundary conditions

t0 = 0.0 # starting point
tf = 10.0 # ending point
N = 1000 # number of points between a and b
dt = (tf-t0)/N

r = np.array([100,0],float) # initial condition

tpoints = np.arange(t0, tf, dt)
ypoints = []
vpoints = []

for t in tpoints:
    ypoints.append(r[0])
    vpoints.append(r[1])
    k1 = dt*f(r,t)
    k2 = dt*f(r+0.5*k1,t+0.5*dt)
    k3 = dt*f(r+0.5*k2,t+0.5*dt)
    k4 = dt*f(r+k3, t+dt)
    r = r + (k1+2*k2+2*k3+k4)/6


# In[6]:


fig0, ax0 = plt.subplots()
ax0.plot(tpoints, ypoints)


# In[ ]:





# In[ ]:





# In[ ]:




