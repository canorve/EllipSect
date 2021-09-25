
from ellipsect.lib.libs import *

from ellipsect import *



def  Interpol(X,Y,X2):

    tck = interpolate.splrep(X, Y)
    return interpolate.splev(X2, tck)

#sectors/num.py
def GetK(n):
    "Solve the Sersic function to get the dependence of K over Sersic index"

    k = gammaincinv(2*n,0.5)

    return (k)


#sectors/num.py
# Deprecated. Do not used  it
def GetKAprox(n):
    "Aproximation to solve the dependence of K on the Sersic index"


    K = 2 * n - 1/3 + 4/(405*n) + 46 / (25515*n**2) + 131 / (1148175 * n**3) - 2194697 / (30690717750*n**4)


    return (K)

def Re90(rad,n):
    "Returns the radius containing the 90% of light"

    #Rad90= galpar.rad * (1.53 + 0.73 * galpar.serind+ 0.07 * galpar.serind**2) 
    x = gammaincinv(2*n,0.9)
    k = gammaincinv(2*n,0.5)

    r9re=(x/k)**n

    rad90 = rad*r9re


    return (rad90)


def RadGamma(rb,alpha,beta,gamma):

    rg = rb * ((1/2 - gamma)/(beta - 1/2))**(1/alpha)

    return rg

