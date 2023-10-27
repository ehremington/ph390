#!/usr/bin/env python
# coding: utf-8

# In[5]:


from scripts import cRK4


# # UGGGGHHHHHHHHHHHHHHHHHHHH

# In[ ]:


def f(r,t):
    x = r[0]
    y = r[1]
    a =
    b = 
    g =
    d =
       
    fx = a*x-b*x*y
    fy = 
    return(np.array([fx,fy],float))

def cRK4(f, tf, x0, v0, t0=0, dt=2**-5): 

    r = np.array([x0,v0],float) # initial condition

    tpoints = np.arange(t0, tf, dt)
    xpoints = []
    vpoints = []

    for t in tpoints:
        xpoints.append(r[0])
        vpoints.append(r[1])
        k1 = dt*f(r,t)
        k2 = dt*f(r+0.5*k1,t+0.5*dt)
        k3 = dt*f(r+0.5*k2,t+0.5*dt)
        k4 = dt*f(r+k3, t+dt)
        r = r + (k1+2*k2+2*k3+k4)/6
        
    return(tpoints, xpoints, vpoints)a

