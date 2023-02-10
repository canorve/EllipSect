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




def GetExpTime(Image,imgidx,flagidx,num,flagnum):
    "Get exposition time from the image"

    hdu = fits.open(Image)
    if flagidx:
        if flagnum:
            exptime = hdu[imgidx,num].header.get("EXPTIME",1) # return 1 if not found
        else:    
            exptime = hdu[imgidx].header.get("EXPTIME",1) # return 1 if not found
    else:
        exptime = hdu[0].header.get("EXPTIME",1) # return 1 if not found

    hdu.close()
    return float(exptime)

#io/fits.py
def GetFits(Image, Imageout, xlo, xhi, ylo, yhi):
    "Get a piece from the image"

    if os.path.isfile(Imageout):
        print("{} deleted; a new one is created \n".format(Imageout))
        runcmd = "rm {}".format(Imageout)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)

    hdu = fits.open(Image)
    #dat = hdu[0].data[ylo - 1: yhi + 1, xlo - 1: xhi + 1]
    dat = hdu[0].data[ylo - 1: yhi, xlo - 1: xhi]
    hdu[0].data = dat
    try:
        hdu.writeto(Imageout, overwrite=True)
    except: 
        hdutemp=fits.PrimaryHDU(data=dat)
        hdutemp.writeto(Imageout, overwrite=True)

    hdu.close()

#io/fits.py
def MakeImage(newfits, sizex, sizey):
    "create a new blank Image"

    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))

        runcmd = "rm {}".format(newfits)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex))
    hdu.writeto(newfits, overwrite=True)

    return True

#io/fits.py
def GetAxis(Image,imgidx,flagidx,num,flagnum):
    # k Check
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    if flagidx:
        if flagnum:
            ncol = hdu[imgidx,num].header.get("NAXIS1",2000) # return 2000 if not found
            nrow = hdu[imgidx,num].header.get("NAXIS2",2000) # return 2000 if not found
        else:    
            ncol = hdu[imgidx].header.get("NAXIS1",2000) # return 2000 if not found
            nrow = hdu[imgidx].header.get("NAXIS2",2000) # return 2000 if not found

    else:
        ncol = hdu[0].header.get("NAXIS1",2000) # return 2000 if not found
        nrow = hdu[0].header.get("NAXIS2",2000) # return 2000 if not found


    hdu.close()

    return ncol, nrow

################################################
################################################
################################################


