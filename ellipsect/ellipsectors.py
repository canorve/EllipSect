#! /usr/bin/env python3



from ellipsect.lib.libs import *


from ellipsect.inout.read import InputSys
from ellipsect.sectors.sect import SectorsGalfit


def run():

    if (len(sys.argv[1:]) == 0):
        print ("EllipSect: an analysis tool for GALFIT output ")
        print ('Missing arguments')
        print ("Usage:\n %s [GALFITOutputFile] [-options] " % (sys.argv[0]))
        print ("use help to display more information about 'options' arguments: ")
        print ("%s -help " % (sys.argv[0]))

        sys.exit()

    # read user's input 
    params = InputSys(sys.argv)

    # full program:

    #photapi stores all the variables computed by EllipSect

    photapi = SectorsGalfit(params)


    #print("AIC: ",photapi.AICrit)
    #print("Bulge to Total: ",photapi.BulgeToTotal)

    return True


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



