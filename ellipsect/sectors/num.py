
from ellipsect.lib.libs import *

from ellipsect import *




def  Interpol(X,Y,X2):

    tck = interpolate.splrep(X, Y)
    return interpolate.splev(X2, tck)

#sectors/num.py
def GetK(n):
    "Solve the Sersic function to get the dependence of K over Sersic index"

## solve the Sersic equation
# to get the dependence of K over
# Sersic index

    count = 1

    #limits
    lima=0
    limb=100

#fx is the function to solve
    fxa = fx(n,lima)
    fxb = fx(n,limb)

    resk= (lima + limb)/2

    fxres=fx(n,resk)


    if(fxa * fxb < 0):

        while(np.abs(fxres) > 0.00000001):

            if(fxa * fxres > 0):
                lima=resk
            elif(fxa * fxres < 0):
                limb=resk
            elif(fxres==0):
                break
            resk= (lima + limb)/2
            fxres=fx(n,resk)

            count+=1

            if (count >= 10000):
                break

    else:
        print("no solution in the range: ({},{})\n".format(lima,limb))

    return (resk)


#sectors/num.py
def fx(n,k):
    "function to solve to get the relation between Sersic index and K"


    func = np.exp(scipy.special.gammaln(2*n)) - 2 * np.exp(scipy.special.gammaln(2*n)) * scipy.special.gammainc(2*n,k)


    return(func)


#sectors/num.py
def GetKAprox(n):
    "Aproximation to solve the dependence of K on the Sersic index"


    K = 2 * n - 1/3 + 4/(405*n) + 46 / (25515*n**2) + 131 / (1148175 * n**3) - 2194697 / (30690717750*n**4)


    return (K)


