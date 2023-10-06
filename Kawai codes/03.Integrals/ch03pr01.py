# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*  Example  3.1                                                          *
%*  filename: ch03pr01.py                                                 *
%*  program listing number: 3.1                                           *
%*                                                                        *
%*  This program integrate sin(x) from x=0 to x=pi using rectangular,     *
%*  trapezoidal and simpson methods.  Absolute errors are plotted.        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/14/2019.                                    *
%**************************************************************************
"""

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

opt=input("Use SciPy [y/n] ")
if opt=='y':
    print("SciPy integrate will be used.\n")
else:
    print("SciPy will not be used.\n")
    
# Set the lower and upper bound of the integration
a=0.
b=np.pi/2.
# Header of the output
print("{0:^75}".format('Absolute error in various numrical integration'))
print("{0:^6} {1:^23} {2:^24} {3:^24} \n"
      .format('N','Rectangular','Trapezoidal','Simpson'))

kmax=10
h=np.zeros(kmax+1)
err_rect=np.zeros(kmax+1)
err_trap=np.zeros(kmax+1)
err_simp=np.zeros(kmax+1)

for k in range(0,kmax):
    N=2**(k+1)
    h[k]=(b-a)/N

    x = a + np.linspace(a,b,N+1)
    f = np.sin(x)
    
    rect=f[0:N].sum()*h[k]

    if opt=='y':
        trap=integrate.trapz(f,x)
        simp=integrate.simps(f,x)
    else:      
        trap=f[1:N].sum()*h[k]+(f[0]+f[N])*h[k]/2.    
        simp=(2.0*f[0:N-1:2].sum()+4.0*f[1:N:2].sum()-f[0]+f[N])*h[k]/3. 
    
    err_rect[k]=abs(1.-rect)
    err_trap[k]=abs(1.-trap)
    err_simp[k]=abs(1.-simp)

    print("{0:5d} {1:24.16e} {2:24.16e} {3:24.16e}"
          .format(N,err_rect[k],err_trap[k],err_simp[k]))
    
    del x
    del f
    
# Plot data
h2=h**2
h3=h**3
h4=h**4

plt.ioff()
plt.figure(figsize=(6,5))
plt.loglog(h,err_rect, 'og', label='rectangular')
plt.loglog(h,err_trap, 'ob', label='trapezoidal')
plt.loglog(h,err_simp, 'or', label='simpson')
plt.loglog(h,h,'--g',label='$h$')
plt.loglog(h,h2,'--b',label='$h^2$')
plt.loglog(h,h4,'--r',label='$h^4$')
plt.legend(loc=4)
plt.xlabel('h')
plt.ylabel('Integral')
plt.show()