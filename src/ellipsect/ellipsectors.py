#! /usr/bin/env python3

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

from ellipsect.inout.prt import printWelcome 

from ellipsect.inout.read import ArgParsing
from ellipsect.sectors.sect import SectorsGalfit


def run():


    printWelcome() # print version and description 


    
    # read user's input 
    args = ArgParsing(sys.argv[1:])

    # full program:
    photapi = SectorsGalfit(args)


    #photapi stores all the variables computed by EllipSect

    #print("AIC: ",photapi.AICrit)
    #print("Bulge to Total: ",photapi.BulgeToTotal)



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


