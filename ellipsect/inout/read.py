
from ellipsect.lib.libs import *

from ellipsect.inout.gfits import GetExpTime
from ellipsect.inout.gfits import GetAxis
from ellipsect.inout.gfits import GetFits



import argparse


def ArgParsing():
    ''' Read user's input '''


    parser = InitParsing()

    args = parser.parse_args()

    #args = parser.parse_args('') # get default parameters


    #params = PassArgs(args)




    return args


def InitParsing():
    '''initial input parsing'''


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
    parser.add_argument("-snr","--snr", action="store_true", help="Creates Signal to Noise image")

    parser.add_argument("-k","--keep", action="store_true", help="use existing file to compute subcomponents")

 
    #options with arguments
    parser.add_argument("-q","--axisrat", type=float, help="galaxy axis ratio ",default=1)
    parser.add_argument("-pa","--posangle", type=float, help="position angle (same as GALFIT)  ",default=0)
    parser.add_argument("-rx","--ranx",nargs=2, type=float, help="provide a range for x-axis: xmin-xmax ")
    parser.add_argument("-ry","--rany", nargs=2,type=float, help="provide a range for y-axis: ymin-ymax  ")
    parser.add_argument("-dpi","--dotsinch", type=int, help="dots per inch used for images files ")
    parser.add_argument("-ml","--minlevel", type=float, help="parameter given directly to sectors_photometry. ")
    parser.add_argument("-sc","--sectors", type=int, help="parameter given directly to sectors_photometry. Divide elipse in 'sectors'")
    parser.add_argument("-ob","--object",  help="used for 'phot' to search in NED")
    parser.add_argument("-f","--filter",  help="used for 'phot' to indicate band for NED  ")
    parser.add_argument("-dm","--distmod", type=float, help="Introduce Distance Modulus  ")
    parser.add_argument("-mc","--magcor", type=float, help="Introduce Galactic Extinction ")
    parser.add_argument("-sk","--scalekpc", type=float, help="Introduce equivalence of ''/kiloparsec ")
    parser.add_argument("-sd","--sbdim", type=float, help="surface brightness dimming ")
    parser.add_argument("-md","--model",  help="User can introduce its own image model. ")
    parser.add_argument("-sky","--sky", type=float, help="User can introduce his/her own sky value. ")
    parser.add_argument("-ned","--ned",  help="user can introduce his/her own ned xml file")
    parser.add_argument("-si","--skyinit", type=float, help="for randsky, it creates a mask for the main target using this radio. For gradsky it is where the program starts to compute the gradient ")
    parser.add_argument("-srm","--skyradmax", type=float, help="for randsky only, maximum radius from main target where randbox can be selected")
    parser.add_argument("-skb","--skybox", type=int, help="pixel size of the box for randsky. Default = 20",default=20)
    parser.add_argument("-skn","--skynum", type=int, help="Number of boxes used in randsky. Default = 20",default=20)
    parser.add_argument("-skw","--skywidth", type=int, help="width of the ring for gradsky. Default = 20 ",default=20)
    parser.add_argument("-distm","--distmax", type=float, help="maximum distance among model centers to be considered to be part of the same galaxy ")
    parser.add_argument("-fw","--fwhm", type=float, help="It is used to compute Area_psf for BICres. Default = 2 pixels ",default=2)


    return parser







#io/read.py
def ReadGALFITout(inputf,galpar,distmax):

    flagfirst = True

    maskimage = ""
    #    skylevel=0

    GalfitFile = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    while index < len(lines):

        #================================================================================
        # IMAGE and GALFIT CONTROL PARAMETERS
        #A) tempfits/A2399-3-2.fits      # Input data image (FITS file)
        #B) A2399-215615.96-m073822.7-337-out.fits      # Output data image block
        #C) tempfits/none-3-2.fits      # Sigma image name (made from data if blank or "none")
        #D) psfs/PSF-1309-721.fits          # Input PSF image and (optional) diffusion kernel
        #E) 1                   # PSF fine sampling factor relative to data
        #F) mask-337            # Bad pixel mask (FITS image or ASCII coord list)
        #G) constraints         # File with parameter constraints (ASCII file)
        #H) 129  809  265  809  # Image region to fit (xmin xmax ymin ymax)
        #I) 60     60           # Size of the convolution box (x y)
        #J) 21.672              # Magnitude photometric zeropoint
        #K) 0.680  0.680        # Plate scale (dx dy)   [arcsec per pixel]
        #O) regular             # Display type (regular, curses, both)
        #P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

        line = lines[index]
        (tmp) = line.split()
        if tmp[0] == "A)":     # input image
            galpar.inputimage=tmp[1]

        if tmp[0] == "B)":     # out image
            galpar.outimage=tmp[1]

        if tmp[0] == "F)":     # mask image
            try:
                galpar.maskimage=tmp[1]
            except IndexError:
                galpar.maskimage="None"


        if tmp[0] == "H)":     # region fit box
            galpar.xmin=int(tmp[1])
            galpar.xmax=int(tmp[2])
            galpar.ymin=int(tmp[3])
            galpar.ymax=int(tmp[4])

        if tmp[0] == "J)":     # mgzpt
            galpar.mgzpt=float(tmp[1])

        if tmp[0] == "K)":     # plate scale
            galpar.scale=float(tmp[1])


        # first component
        if tmp[0] == "0)" and flagfirst == True:     # plate scale

            flagfirst=False

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    galpar.xc=float(tmp[1])
                    galpar.yc=float(tmp[2])

                if tmp[0] == "4)":    # Effective radius
                    galpar.rad=float(tmp[1])

                if tmp[0] == "5)":    # Sersic index 
                    galpar.serind=float(tmp[1])


                if tmp[0] == "9)":    # axis ratio
                    galpar.q=float(tmp[1])

                if tmp[0] == "10)": # position angle
                    galpar.ang=float(tmp[1])

        # sersic component
        if tmp[0] == "0)" and tmp[1] == "sersic":     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcser=float(tmp[1])
                    ycser=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((galpar.xc-xcser)**2+(galpar.yc-ycser)**2)


                if tmp[0] == "4)" and (dist < distmax):    # Effective radius
                    galpar.rad=float(tmp[1])

                if tmp[0] == "5)" and (dist < distmax):    # Sersic index 
                    galpar.serind=float(tmp[1])



                if tmp[0] == "9)" and (dist < distmax ):    # axis ratio
                    galpar.q=float(tmp[1])

                if tmp[0] == "10)" and (dist < distmax): # position angle
                    galpar.ang=float(tmp[1])

        # second component exponential model
        if tmp[0] == "0)" and tmp[1] == "expdisk" :     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcexp=float(tmp[1])
                    ycexp=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((galpar.xc-xcexp)**2+(galpar.yc-ycexp)**2)


                if tmp[0] == "4)" and (dist < distmax):    # Effective radius
                    galpar.rad=float(tmp[1])
                    galpar.rad=1.678*galpar.rad

                #if tmp[0] == "5)" and (dist < distmax):    # Sersic index 
                    galpar.serind=1



                if (tmp[0] == "9)") and (dist < distmax):    # axis ratio
                    galpar.q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < distmax): # position angle
                    galpar.ang=float(tmp[1])


        if tmp[0] == "0)" and tmp[1] == "gaussian":

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcgauss=float(tmp[1])
                    ycgauss=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((galpar.xc-xcgauss)**2+(galpar.yc-ycgauss)**2)


                if tmp[0] == "4)" and (dist < distmax):    # Effective radius
                    galpar.rad=float(tmp[1])
                    galpar.rad=0.8325546*galpar.rad

                #if tmp[0] == "5)":    # Sersic index 
                    galpar.serind=0.5



                if (tmp[0] == "9)") and (dist < distmax):    # axis ratio
                    galpar.q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < distmax): # position angle
                    galpar.ang=float(tmp[1])


        if tmp[0] == "0)" and tmp[1] == "sky" :     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":    # axis ratio
                    galpar.skylevel=float(tmp[1])

        index += 1

    GalfitFile.close()

    ###################
    chars = set('[]') 
    numbers=set('1234567890')
    if any((c in chars) for c in galpar.inputimage): 
        print("Ext Found") 
        galpar.flagidx=True
        (filename,imgidxc) = galpar.inputimage.split("[")
        (imgidx,trash)=imgidxc.split("]")

        if any((n in numbers) for n in imgidx):
            galpar.flagnum=True
            (imgidx,num)=imgidx.split(",")
            num=int(num)

        galpar.inputimage=filename
        galpar.imgidx=imgidx
        galpar.num=num

    #    else: 
    #       print("Not found") 
    
    ####################

    errmsg="file {} does not exist".format(galpar.inputimage)
    assert os.path.isfile(galpar.inputimage), errmsg

    galpar.exptime=GetExpTime(galpar.inputimage,galpar.imgidx,galpar.flagidx,galpar.num,galpar.flagnum)


    #errmsg="xc and yc are unknown "
    #assert ("xc" in locals()) and ("yc" in locals())  , errmsg

    print("center is at xc, yc = ",galpar.xc,galpar.yc)

    galpar.inxc=galpar.xc
    galpar.inyc=galpar.yc

    # correcting coordinates
    galpar.xc=galpar.xc-galpar.xmin+1
    galpar.yc=galpar.yc-galpar.ymin+1


    flagaxis=False
    if (os.path.isfile(galpar.maskimage)):
        col1,row1=GetAxis(galpar.maskimage,"none",False,1,False)
    else:
        col1=row1=0
    
    col2,row2=GetAxis(galpar.inputimage,galpar.imgidx,galpar.flagidx,galpar.num,galpar.flagnum)

    if (col1 != col2 or row1 != row2):
        flagaxis=True


    if (os.path.isfile(galpar.maskimage)) and (flagaxis==False):
        mime=mimetypes.guess_type(galpar.maskimage)

        flagbm = not(mime[0] == "text/plain")

        errmsg="Sorry the mask file: {}  must be binary, not ASCII ".format(maskimage)
        assert flagbm is True, errmsg

        GetFits(galpar.maskimage, galpar.tempmask, galpar.xmin, galpar.xmax, galpar.ymin, galpar.ymax)

    else:
        if flagaxis == True:
            errmsg="Mask and image do not have the same shape. Ignoring mask"
            print(errmsg)
        else:
            errmsg="Unable to find Mask file"
            print(errmsg)


        
        #creates empty mask
        hdu = fits.open(galpar.outimage)
        tempimg = (hdu[1].data.copy()).astype(float)
        hdu.close()

        yaxis=tempimg.shape[0]
        xaxis=tempimg.shape[1]

        data = np.zeros((yaxis , xaxis ), dtype=np.float64)
        hdu = fits.PrimaryHDU(data=data)
        hdu.writeto(galpar.tempmask, overwrite=True)




    #return xc,yc,q,pa,skylevel,scale,outimage,mgzpt,exptime,mask

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



