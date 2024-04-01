import numpy as np
from scipy import integrate
from numpy import exp
from math import sqrt
import math 

def f(t,s):
    a = 4*math.pi**2
    m = ((4-t**2)**0.5)*((4-s**2)**0.5)
    l = 4-t*s
    f1 = 1
    f2 = 1
    return (a**-1)*f1*f2*np.log((l+m)/(l-m))

v, err = integrate.dblquad(f,-2,2,lambda g : -2,lambda h : 2)
