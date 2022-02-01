#! /usr/bin/env python3



from ellipsect.lib.libs import *


from ellipsect.inout.read import ArgParsing
from ellipsect.sectors.sect import SectorsGalfit


def run():

    # read user's input 
    args = ArgParsing()

    # full program:

    #photapi stores all the variables computed by EllipSect

    photapi = SectorsGalfit(args)


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



