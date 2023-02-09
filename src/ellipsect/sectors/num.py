#    EllipSect: An analysis tool for GALFIT output 
#    Copyright (C) 2022  Christopher Añorve 

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


from ellipsect.lib.libs import *

from ellipsect import *



def  Interpol(X,Y,X2):

    tck = interpolate.splrep(X, Y)

    flag = False
    if(np.isnan(np.sum(tck[1]))):

        m=len(X)
        s=m-np.sqrt(2*m)
        tck = interpolate.splrep(X, Y,s=s)
        flag=True

    Y2 = interpolate.splev(X2, tck)

    return Y2,flag

#sectors/num.py
def GetK(n):
    "Solve the Sersic function to get the dependence of K over Sersic index"

    k = sc.gammaincinv(2*n,0.5)

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
    x = sc.gammaincinv(2*n,0.9)
    k = sc.gammaincinv(2*n,0.5)

    r9re=(x/k)**n

    rad90 = rad*r9re


    return (rad90)


def RadGamma(rb,alpha,beta,gamma):

    rg = rb * ((1/2 - gamma)/(beta - 1/2))**(1/alpha)

    return rg

def ReFrac(rad,n,frac):
    "Returns the radius containing the fraction of light"

    x = sc.gammaincinv(2*n,frac)
    k = sc.gammaincinv(2*n,0.5)

    rfre=(x/k)**n

    radfrac = rad*rfre

    return (radfrac)



def KronRadius(rad,re,n):
    "return the Kron Radius"

    k = sc.gammaincinv(2*n,0.5)

    x = k*(rad/re)**(1/n)

    KronRad = ((re/k**n))*((sc.gammainc(3*n,x)*sc.gamma(3*n))/(sc.gammainc(2*n,x)*sc.gamma(2*n)))

    return KronRad


def PetrosianIndex(rad,re,n,C=0):
    "returns the inverse of the Petrosian Index"

    k = sc.gammaincinv(2*n,0.5)

    x = k*(rad/re)**(1/n)


    nu = (2*n*(sc.gammainc(2*n,x)*sc.gamma(2*n)))/((np.exp(-x))*x**(2*n))

    invnu  = 1/nu

    return invnu - C


def solvePet(n,c=0.2):
    "return the Petrosian radius (in Re units). It uses Bisection"


    a = .1
    b = 30

    if(isinstance(n,np.ndarray)):
        Rp =[]

        for idx,item in enumerate(n):
            r= bisect(PetrosianIndex,a,b,args=(1,item,c))
            Rp.append(r)

    else:

        Rp=bisect(PetrosianIndex,a,b,args=(1,n,c))

    return Rp



