#!/usr/bin/python

"""
This is a file that contains functions that we will use over and over again.

"""

def cRK4(f, tf, x0, v0, t0=0, dt=2**-5): 
    from numpy import array, arange
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
        
    return(tpoints, xpoints, vpoints)

if __name__ == "__main__":
    print('you didn\'t mean to do this')