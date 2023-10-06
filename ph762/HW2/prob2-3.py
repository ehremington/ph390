import calculus as calc
import numpy as np
def f(x):
    return x*np.cos(x)

# here h is the exponent on a base of 2
h = 20

#test at a variety of steps sizes
for i in range(h):
    steps = 2**i+1

    FT = calc.IntTrapFromFunction(f,0.0,np.pi,steps)

    FS = calc.IntSimpFromFunction(f,0.0,np.pi,steps)

    print('%i      '%steps + '%.8f         '%FT+'%.8f        '%FS)