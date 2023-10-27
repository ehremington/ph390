#!/usr/bin/env python
# coding: utf-8

# # Day 29

# In[6]:


def f(r,t):
    x = r[0]
    y = r[1]
    z = r[2]
    s = 10
    t = 28
    b = 8/3
    fx = s*(y-x)
    fy = t*x-y-x*z
    fz = x*y-b*z
    return(np.array([fx, fy, fz],float))

def cRK4(f, tf, x0, y0, z0, t0=0, dt=2**-5): 

    r = np.array([x0,y0,z0],float) # initial condition

    tpoints = np.arange(t0, tf, dt)
    xpoints = []
    ypoints = []
    zpoints = []

    for t in tpoints:
        xpoints.append(r[0])
        ypoints.append(r[1])
        zpoints.append(r[2])
        k1 = dt*f(r,t)
        k2 = dt*f(r+0.5*k1,t+0.5*dt)
        k3 = dt*f(r+0.5*k2,t+0.5*dt)
        k4 = dt*f(r+k3, t+dt)
        r = r + (k1+2*k2+2*k3+k4)/6
        
    return(tpoints, xpoints, ypoints, zpoints)


# In[7]:


t, x, y, z = cRK4(f, 50, 0, 1, 0)


# In[8]:


fig0, ax0 = plt.subplots()
ax0.plot(t,y)


# In[10]:


fig1, ax1 = plt.subplots()
ax1.plot(x, z)


# In[ ]:


def projectile(f, tf, x0, y0, z0, t0=0, dt=2**-5): 

    r = np.array([x0,y0,z0],float) # initial condition

    tpoints = np.arange(t0, tf, dt)
    xpoints = []
    ypoints = []
    zpoints = []

    for t in tpoints:
        xpoints.append(r[0])
        ypoints.append(r[1])
        zpoints.append(r[2])
        k1 = dt*f(r,t)
        k2 = dt*f(r+0.5*k1,t+0.5*dt)
        k3 = dt*f(r+0.5*k2,t+0.5*dt)
        k4 = dt*f(r+k3, t+dt)
        r = r + (k1+2*k2+2*k3+k4)/6
        
    return(tpoints, xpoints, ypoints, zpoints)


# In[ ]:


def f(r,t):
    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]
    R = 
    m = 
    g = 9.8
    rho = 
    C = 
    fx = 
    fvx = 
    fy = 
    fvy = 
    
    return(np.array([fx, fy, fz],float))

