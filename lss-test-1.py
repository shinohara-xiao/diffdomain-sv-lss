import numpy as np
import math
from scipy import stats
from scipy import integrate

def g(x):
    #return x**2
    e = np.exp(x)
    return e

def Ff(x):
    out = np.sqrt(4 - x**2) * g(x) / 2 / math.pi
    return out
    
def lss_formula(Mat):
    nd = Mat.shape[0]
    # normalize it to a Wigner matrix
    Mat = Mat / np.sqrt(nd)
    try:
        w, v = np.linalg.eigh(Mat)
        #lambdan = np.amax(np.abs(w))
        print('max lambdan:', np.amax(np.abs(w)))
        
        # calculate the Gn(f)
        GnF = integrate.quad(Ff,-2,2)[0]
        GnF = GnF * nd   
        
        print('1')
        
        a = np.sum(g(w))
        Gnf = a - GnF
        
        print('3')
        
        lambdan = (Gnf - 1)/np.sqrt(4)
        
        
        if (Gnf-1) >= 0:
            p = 1 - stats.norm.cdf(lambdan)
            
        else:
            p = stats.norm.cdf(lambdan)
            
        p = 2 * p
        
        # normalize by sqrt(n)
        #lambdan =  np.power(nd, 2.0/3) * (lambdan -2)

        #twdist = tw.TracyWidom(beta=1)
        #p = 1 - twdist.cdf(lambdan)
    except:
        nd = np.nan
        lambdan = np.nan
        p = np.nan
        
    return nd, lambdan, p 

matrix = np.mat(np.random.rand(20,20))

test = lss_formula(matrix)
test

# matrix = np.mat(np.random.randn(20,20))
def g(x):
    return np.exp(x)

matrix = np.mat(np.random.rand(20,20))
w, v = np.linalg.eigh(matrix)
g(w)

def Ff(x):
    out = np.sqrt(4-x**2)*g(x)/2/math.pi
    return out
Ff(1)

nd = matrix.shape[0]

def Ff(x):
    out = np.sqrt(4 - x**2) * g(x) / 2 / math.pi
    return out

GnF = integrate.quad(Ff,-2,2)[0]
GnF = GnF * nd

a = np.sum(g(w))
Gnf = a - GnF
lambdan = (Gnf - 0)/np.sqrt(2)

if lambdan >= 0:
    p = 1 - stats.norm.cdf(lambdan)
else:
    p = stats.norm.cdf(lambdan)

p = 2 * p

p
