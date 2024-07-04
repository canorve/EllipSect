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

from ellipsect.inout.gfits import GetExpTime
from ellipsect.inout.gfits import GetAxis
from ellipsect.inout.gfits import GetFits


from ellipsect.inout.galfit import SelectGal 


from astropy.wcs import WCS

import argparse


def ArgParsing(args):
    ''' Read user's arguments input '''


    parser = InitParsing()

    parser_args = parser.parse_args(args)



    return parser_args


def InitParsing():
    '''initial argument parsing'''

    # every new parameter added here must be
    # also added in function PassArgs() in ellipsect/sectors/sect.py


    parser = argparse.ArgumentParser(description="EllipSect creates surface brightness profiles (and other photometric variables) from GALFIT output ")

    # required arguments
    parser.add_argument("GalFile",default="galfit.01",help="GALFIT output file: galfit.XX ")

    #options without arguments
    parser.add_argument("-lx","--logx", action="store_true", help="turn the X-axis to logarithm ")
    parser.add_argument("-cp","--comp", action="store_true", help="add individual model components to the plot")
    parser.add_argument("-px","--pix", action="store_true", help="turn the top x-axis in pixels ")
    parser.add_argument("-g","--grid", action="store_true", help="display a grid in the plot ")
    parser.add_argument("-sb","--sbout", action="store_true", help="creates output file containing the surface brightness profiles")
    parser.add_argument("-np","--noplot", action="store_true", help="avoids pop up windows and only creates images files")
    parser.add_argument("-ph","--phot", action="store_true", help="Compute photometry. Check the created output file")
    parser.add_argument("-ci","--checkimg", action="store_true", help="save the images used for sectors_photometry in individual files")
    parser.add_argument("-nn","--noned", action="store_true", help="it avoids to connect to NED")
    parser.add_argument("-gsky","--gradsky", action="store_true", help="computes sky using the gradient method")
    parser.add_argument("-rsky","--randsky", action="store_true", help="computes sky using random boxes")
    parser.add_argument("-snr","--snr", action="store_true", help="Creates Signal to Noise image if phot is activated")
    #parser.add_argument("-chi","--chisqr", action="store_true", help="Creates Chi-square image if phot is activated")


    parser.add_argument("-re","--effrad", action="store_true", help="Draw a vertical line indicating the effective radius")

    parser.add_argument("-r90","--rad90", action="store_true", help="Draw a vertical line indicating the 90%% of total light")

    parser.add_argument("-r95","--rad95", action="store_true", help="Draw a vertical line indicating the 95%% of total light")


    parser.add_argument("-rc","--remcomp", action="store_true", help="remove file of subcomponents")
    parser.add_argument("-gx","--galax", action="store_true", help="only the galaxy surface brightness is shown in the plot")


    linehelp="includes all the sky pixels for gradsky and randsky. "  \
            + "Default mode removes top %%80 and bottom %%20 sky pixels "
    parser.add_argument("-asp","--allskypx", action="store_true", 
                        help=linehelp)



    parser.add_argument("-t","--title", action="store_true", help="put title on plots")

    #options with arguments
    #parser.add_argument("-cn","--center",nargs=2,action="store", type=float, help="galaxy's center ")

    parser.add_argument("-q","--axisrat", type=float, help="galaxy axis ratio ")
    parser.add_argument("-pa","--posangle", type=float, help="position angle (same as GALFIT)  ",default=0)
    parser.add_argument("-rx","--ranx",nargs=2, type=float, help="provide a range for x-axis: xmin - xmax ")
    parser.add_argument("-ry","--rany", nargs=2,type=float, help="provide a range for y-axis: ymin - ymax  ")
    parser.add_argument("-dpi","--dotsinch", type=int, help="dots per inch used for images files ")
    parser.add_argument("-ml","--minlevel", type=float, help="parameter given to sectors_photometry. ")
    parser.add_argument("-sc","--sectors", type=int, help="parameter given to sectors_photometry. It divides the ellipse in 'sectors'")
    parser.add_argument("-ob","--object",  help="used for 'phot' to search in NED")
    parser.add_argument("-f","--filter",  help="used for 'phot' to indicate band for NED  ")
    parser.add_argument("-dm","--distmod", type=float, help="Introduce Distance Modulus  ")
    parser.add_argument("-mc","--magcor", type=float, help="Introduce Galactic Extinction ")
    parser.add_argument("-sk","--scalekpc", type=float, help="Introduce equivalence of ''/kiloparsec ")
    parser.add_argument("-sd","--sbdim", type=float, help="surface brightness dimming ")
    parser.add_argument("-md","--model",  help="User can introduce its own image model. ")
    parser.add_argument("-sky","--sky", type=float, help="Program will use this sky values instead of the GALFIT")
    parser.add_argument("-ned","--ned",  help="user can introduce his/her own ned xml file")
    parser.add_argument("-ri","--radinit", type=float, help="for randsky, it creates a mask for the main target using this radio. For gradsky it is where the program starts to compute the gradient ")
    parser.add_argument("-srm","--skyradmax", type=float, help="for randsky only, maximum radius from main target where randbox can be selected")
    parser.add_argument("-skb","--skybox", type=int, help="pixel size of the box for randsky. Default = 20",default=20)
    parser.add_argument("-skn","--skynum", type=int, help="Number of boxes used in randsky. Default = 20",default=20)
    parser.add_argument("-skw","--skywidth", type=int, help="width of the ring for gradsky. Default = 20 ",default=20)
    parser.add_argument("-distm","--distmax", type=float, help="maximum distance among model centers to be considered to be part of the same galaxy ")
    parser.add_argument("-fw","--fwhm", type=float, help="It is used to compute Area_psf for BICres. It also draws a vertical line at the given value. Default = 2 pixels ")


    #parser.add_argument("-fr","--frac", type=float, help="fraction of the minimum model count to be used as vmin of imshow for the cube image ",default=1)

    #parser.add_argument("-frm","--fracmax", type=float, help="fraction of the maximum model count to be used as vmax of imshow for the cube image ",default=1)


    parser.add_argument("-br","--brightness", type=float, 
                        help="brightness of the image. Only for galaxy and model. Default = 0. Preferible range goes from -1 to 1", default=0)
    parser.add_argument("-co","--contrast", type=float, 
                        help="contrast of the image. Only for galaxy and model. Default = 1. Preferible range goes from 0 to 1",default=1)



    parser.add_argument("-cm","--cmap", type=str, help="cmap to be used for the cube image ",default="viridis")


    parser.add_argument("-nc","--numcomp", type=int, help="component number to" 
                        + "select for galaxy center. Default = first component. The "
                        + "component order follows as it is shown in galfit file ",default=1)


    parser.add_argument("-ae","--aext", type=float, 
                        help="Surface brightness correction for plots only ", default=0)
 
    parser.add_argument("-hc","--hconst", type=float, 
                        help="hubble constant to download xml file from NED ", default=67.8)
 
    parser.add_argument("-om","--omegam", type=float, 
                        help="omega matter value to download xml file from NED ", default=.308)
 
    parser.add_argument("-ov","--omegav", type=float, 
                        help="omega lambda value to download xml file from NED ", default=0.692)
 

    # every new parameter added here must be
    # also added in function PassArgs() in ellipsect/sectors/sect.py


    return parser





#io/read.py
def ReadGALFITout(ellconf,galhead,galcomps):


    galcomps = SelectGal(galcomps,ellconf.distmax,ellconf.numcomp)
 

    errmsg="file {} does not exist".format(galhead.inputimage)
    assert os.path.isfile(galhead.inputimage), errmsg

    galhead.exptime = GetExpTime(galhead.inputimage,galhead.imgidx,galhead.flagidx,
                                    galhead.num,galhead.flagnum)


    GetEllInfo(ellconf, galcomps) 


    #saving coordinates for the large images: input image, mask img, seg img
    ellconf.inxc = ellconf.xc
    ellconf.inyc = ellconf.yc

    #coordinates for the output images: galaxy, model, residual
    ellconf.xc = ellconf.xc - galhead.xmin + 1
    ellconf.yc = ellconf.yc - galhead.ymin + 1


    flagaxis = False
    if (os.path.isfile(galhead.maskimage)):
        col1, row1 = GetAxis(galhead.maskimage,"none",False,1,False)
    else:
        col1 = row1 =0
    
    col2, row2 = GetAxis(galhead.inputimage, galhead.imgidx, galhead.flagidx,
                        galhead.num, galhead.flagnum)

    if (col1 != col2 or row1 != row2):
        flagaxis = True

    if (os.path.isfile(galhead.maskimage)) and (flagaxis == False):
        mime=mimetypes.guess_type(galhead.maskimage)

        flagbm = not(mime[0] == "text/plain")


        if flagbm is False:

            print("Converting mask ascii to mask FITS ")
            Value = 100  # default value for fits mask
            galhead.maskimage = xy2fits().MakeFits(galhead.inputimage, 
                                                    galhead.maskimage, Value)

        GetFits(galhead.maskimage, galhead.tempmask, galhead.xmin, galhead.xmax, 
                galhead.ymin, galhead.ymax)

    else:
        if flagaxis == True:
            errmsg="Mask and image do not have the same shape. Ignoring mask"
            print(errmsg)
        else:
            errmsg="Unable to find Mask file"
            print(errmsg)

        
        #creates empty mask
        hdu = fits.open(galhead.outimage)
        tempimg = (hdu[1].data.copy()).astype(float)
        hdu.close()

        yaxis = tempimg.shape[0]
        xaxis = tempimg.shape[1]

        data = np.zeros((yaxis, xaxis ), dtype=np.float64)
        hdu = fits.PrimaryHDU(data = data)
        hdu.writeto(galhead.tempmask, overwrite=True)


def GetEllInfo(ellconf,galcomps):
    '''Gets geometry information from the last component''' 


    maskactive = (galcomps.Active == True) 

    if ellconf.flagq == False:
        ellconf.qarg = galcomps.AxRat[maskactive][-1]

    if ellconf.flagpa == False:
        ellconf.parg = galcomps.PosAng[maskactive][-1]

    #for the center, it obtains it from the first component
    #if ellconf.flagcent == False:
    ellconf.xc = galcomps.PosX[maskactive][0]
    ellconf.yc = galcomps.PosY[maskactive][0]



    return True


#inout/read.py
def ReadNComp(inputf,X,Y,galcomps,distmax):
    ## search and count model components

    GalfitFile = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    n=0

    distmin=distmax # minimum distance for among component centers

    while index < len(lines):

        line = lines[index]
        (tmp) = line.split()

        if tmp[0] == "H)":     # region fit box
            xmin=int(tmp[1])
            xmax=int(tmp[2])
            ymin=int(tmp[3])
            ymax=int(tmp[4])

            X=X+xmin-1
            Y=Y+ymin-1

        #init values
        Comps=False
        N=0
        NameComp="none"
        PosX=0
        PosY=0
        Mag=99
        Rad=0
        Exp=0
        Exp2=0
        Exp3=0
        AxRat=1
        PosAng=0
        skip=1
        flagcomp=False  # detect components
        freepar=0
        if tmp[0] == "0)":

            namec=tmp[1] 
            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xc=float(tmp[1])
                    yc=float(tmp[2])

                    dist = np.sqrt((xc-X)**2+(yc-Y)**2)
                    if (dist < distmin and namec != "sky"):
                        n=n+1
                        PosX=xc
                        PosY=yc
                        Comps=True
                        NameComp=namec
                        N= n 
                        flagcomp=True
                        par1=int(tmp[3])
                        par2=int(tmp[4])
                        # count the number of free params
                        freepar=freepar+par1+par2  
                    else:
                        Comps=False

                if tmp[0] == "3)" and flagcomp == True:    # axis ratio
                    Mag=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "4)" and flagcomp == True:    # axis ratio
                    Rad=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "5)" and flagcomp == True:    # axis ratio
                    Exp=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "6)" and flagcomp == True:    # axis ratio
                    Exp2=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "7)" and flagcomp == True:    # axis ratio
                    Exp3=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "9)" and flagcomp == True:    # axis ratio
                    AxRat=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "10)" and flagcomp == True:    # position angle 
                    PosAng=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "z)" and flagcomp == True:    # skip
                    skip=int(tmp[1])
            if (flagcomp == True):

                galcomps.PosX=np.append(galcomps.PosX,PosX)
                galcomps.PosY=np.append(galcomps.PosY,PosY)
                galcomps.Comps=np.append(galcomps.Comps,Comps)
                galcomps.NameComp=np.append(galcomps.NameComp,NameComp)
                galcomps.N=np.append(galcomps.N, N)
                
                galcomps.Mag=np.append(galcomps.Mag,Mag)
                galcomps.Rad=np.append(galcomps.Rad,Rad)
                galcomps.Exp=np.append(galcomps.Exp,Exp)
                galcomps.Exp2=np.append(galcomps.Exp2,Exp2)
                galcomps.Exp3=np.append(galcomps.Exp3,Exp3)
                galcomps.AxRat=np.append(galcomps.AxRat,AxRat)
                galcomps.PosAng=np.append(galcomps.PosAng,PosAng)
                galcomps.skip=np.append(galcomps.skip,skip)
                galcomps.freepar=np.append(galcomps.freepar,freepar)
      
        index += 1

    GalfitFile.close()

    tot=galcomps.Comps.size

    # computed parameters:
    galcomps.Rad50=np.array([0.0]*tot)
    galcomps.SerInd=np.array([0.0]*tot)
    galcomps.Rad50kpc=np.array([0.0]*tot)
    galcomps.Rad50sec=np.array([0.0]*tot)
    galcomps.Rad90=np.array([0.0]*tot)
    galcomps.AbsMag=np.array([99.0]*tot)
    galcomps.Lum=np.array([0.0]*tot)
    galcomps.Flux=np.array([0.0]*tot)
    galcomps.PerLight=np.array([0.0]*tot)
    galcomps.me=np.array([99.0]*tot)
    galcomps.mme=np.array([99.0]*tot)
    galcomps.kser = np.array([0.0]*tot)

    galcomps.KronRad=np.array([0.0]*tot)
    galcomps.PetRad=np.array([0.0]*tot)


    galcomps.N=galcomps.N.astype(int)

    return True


class xy2fits:
    ''' class function that converts the mask ascii file to fits file'''

    def MakeFits(self,ImageFile, AsciiFile, Value):
        
        (tmp)=AsciiFile.split(".")

        namefile=tmp[0]

        maskfits=namefile + ".fits"


        (ncol, nrow) = self.GetAxis(ImageFile)

        self.MakeImage(maskfits, ncol, nrow)

        X, Y = np.genfromtxt(AsciiFile, delimiter="", unpack=True)

        X = X.astype(int)
        Y = Y.astype(int)

        X = X-1
        Y = Y-1


        self.PutPix(X,Y,Value,maskfits)

        print("Ascii -> Fits done ")

        return maskfits

    def GetAxis(self,Image):
        "Get number of rows and columns from the image"

        hdu = fits.open(Image)
        ncol = hdu[0].header["NAXIS1"]  # for hubble images
        nrow = hdu[0].header["NAXIS2"]
        hdu.close()

        return ncol, nrow

    def MakeImage(self,newfits, sizex, sizey):
        "create a new blank Image"

        if os.path.isfile(newfits):
            print("{} deleted; a new one is created ".format(newfits))

            runcmd = "rm {}".format(newfits)
            errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


        hdu = fits.PrimaryHDU()
        hdu.data = np.zeros((sizey, sizex),dtype=np.float64)
        hdu.writeto(newfits, overwrite=True)

        return True



    def PutPix(self,X,Y,Value,ImageFits):

        # original file
        hdu=fits.open(ImageFits)
        Image = hdu[0].data

        ## for some strange reason I have to interchange X and Y
        Image[[Y],[X]]=Value

        hdu[0].data=Image
        hdu.writeto(ImageFits,overwrite=True)
        hdu.close()


def GetWCS(Image):
    "Get World Cordinate System info"
    hdu = fits.open(Image)[0]
    wcs = WCS(hdu.header)

    return wcs 


def prefixNames(ellconf, outimage):
    '''names of the output files based on prefix of galfit output'''

    #note: make one function to save all the names: 
    root_ext = os.path.splitext(outimage)
    ellconf.namefile = root_ext[0]
    
    #adding numbers of galfit file for output
    gal_ext = os.path.splitext(ellconf.galfile)
    numstr = gal_ext[1] 
    numstr = numstr[1:]

    if ellconf.numcomp:
        precomp = "-nc" + str(ellconf.numcomp)
    else:
        precomp = ""




    if numstr.isnumeric():

        # names for the different png
        ellconf.namepng = ellconf.namefile + '-g' + numstr + precomp + "-splot.png"
        ellconf.namesec = ellconf.namefile + '-g' + numstr + precomp + "-gal.png"
        ellconf.namemod = ellconf.namefile + '-g' + numstr + precomp + "-mod.png"
        ellconf.namemul = ellconf.namefile + '-g' + numstr + precomp + "-mplot.png"
        ellconf.namesub = ellconf.namefile + '-g' + numstr + precomp + "-comp.fits"

        ellconf.namesig = ellconf.namefile + '-g' + numstr + precomp + "-sig.fits"


        ellconf.sboutput = ellconf.namefile + '-g' + numstr + precomp + "-sbout"
        ellconf.output = ellconf.namefile + '-g' + numstr + precomp +  "-out.txt"

        ellconf.namened = ellconf.namefile + '-g' + numstr + precomp + "-ned.xml"

        ellconf.namesnr = ellconf.namefile + '-g' + numstr + precomp + "-snr.fits"
        ellconf.namecheck = ellconf.namefile + '-g' + numstr + precomp + "-check.fits"
        ellconf.namering = ellconf.namefile + '-g' + numstr + precomp + "-ring.fits"
        ellconf.nameringmask = ellconf.namefile + '-g' + numstr + precomp + "-ringmask.fits"
        ellconf.namecube = ellconf.namefile + '-g' + numstr + precomp + "-cub.png"
        ellconf.namechi = ellconf.namefile + '-g' + numstr + precomp + "-chi.fits"

    else:

        # names for the different png
        ellconf.namepng = ellconf.namefile + "-splot.png"
        ellconf.namesec = ellconf.namefile + "-gal.png"
        ellconf.namemod = ellconf.namefile + "-mod.png"
        ellconf.namemul = ellconf.namefile + "-mplot.png"
        ellconf.namesub = ellconf.namefile + "-comp.fits"

        ellconf.namesig = ellconf.namefile + "-sig.fits"


        ellconf.sboutput = ellconf.namefile + "-sbout"
        ellconf.output = ellconf.namefile + "-out.txt"

        ellconf.namened = ellconf.namefile + "-ned.xml"

        ellconf.namesnr = ellconf.namefile + "-snr.fits"
        ellconf.namecheck = ellconf.namefile + "-check.fits"
        ellconf.namering = ellconf.namefile + "-ring.fits"
        ellconf.nameringmask = ellconf.namefile + "-ringmask.fits"
        ellconf.namecube = ellconf.namefile + "-cub.png"
        ellconf.namechi = ellconf.namefile + "-chi.fits"




    return True

