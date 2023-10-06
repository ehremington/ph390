#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*  Secion 4.6.                                                           *
%*  filename: ch04pr06.py                                                 *
%*  program listing number: 4.6-1                                         *
%*                                                                        *
%*     Require: findaroot.py                                              *
%*                                                                        *
%*     This program finds trhe energy eigenvalues of a pqrticle           *
%*     in a finite square well potential using the secant root            *
%*     finding methods.                                                   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
"""

import numpy as np
import rootfinding as rf

def f(z):
    global z0
    return z*np.tan(z) - np.sqrt(z0*z0-z*z)

if __name__ == "__main__":
    global z0
    z0=6  #system parameter
    #control parameters
    N=100
    K=np.int(np.ceil(z0/np.pi))  #The number of roots

    # initial bracket
    z1=0.0;
    z2=np.pi/2.0

    for k in range(1,K+1):  
        dz = (z2-z1)/N  # small shit
        z = rf.findaroot(f,z1+dz,z2-dz,10,1.0e-6)
        print("root={0:12.8f}".format(z))
        z1 = z2
        z2 = min([z0,z1 + np.pi])

