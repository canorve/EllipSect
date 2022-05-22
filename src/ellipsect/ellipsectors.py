#! /usr/bin/env python3



from ellipsect.lib.libs import *


from ellipsect.inout.read import ArgParsing
from ellipsect.sectors.sect import SectorsGalfit


def run():

    # read user's input 
    args = ArgParsing(sys.argv[1:])

    # full program:

    photapi = SectorsGalfit(args)


    #photapi stores all the variables computed by EllipSect

    #print("AIC: ",photapi.AICrit)
    #print("Bulge to Total: ",photapi.BulgeToTotal)

    return photapi 


    ##############       #############
    ##############  END  #############
    ##############       #############

    #     ______________________________________________________________________
    #    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
    #   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
    #   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
    #   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
    #   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
    #   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
    #   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/



