
from ellipsect.lib.libs import *

from ellipsect.inout.gfits import GetExpTime
from ellipsect.inout.gfits import GetAxis
from ellipsect.inout.gfits import GetFits

from ellipsect.inout.help import Help
from ellipsect.lib.clas import InputParams

def InputSys(argv):
    ''' Read user's input '''

    OptionHandleList = ['-logx', '-q', '-pa','-comp','-pix','-ranx','-rany','-grid','-dpi','-sbout','-noplot',
        '-minlevel','-sectors','-phot','-object','-filter','-snr','-help','-checkimg','-noned','-distmod','-magcor',
        '-scalekpc','-sbdim','-model','-sky','-keep','-ned','-gradsky','-randsky','-skyRad','-skyRadmax','-skynum','-skybox','-skywidth','-distmax', '-fwhm']


    #class for user's parameters
    params=InputParams()

    options = {}
    for OptionHandle in OptionHandleList:
        options[OptionHandle[1:]] = argv[argv.index(OptionHandle)] if OptionHandle in argv else None
    if options['logx'] != None:
        params.flaglogx=True
        print("X axis is logarithm")
    if options['q'] != None:
        params.flagq=True
    if options['pa'] != None:
        params.flagpa=True
    if options['pix'] != None:
        params.flagpix=True
    if options['ranx'] != None:
        params.flagranx[0]=True
    if options['rany'] != None:
        params.flagrany[0]=True
    if options['grid'] != None:
        params.flagrid=True
    if options['dpi'] != None:
        params.flagdpi=True
    if options['comp'] != None:
        params.flagcomp=True
        print("Plotting subcomponents ")
    if options['noplot'] != None:
        params.flagnoplot=True
        print("images will not be displayed")
    if options['sbout'] != None:
        params.flagsbout=True
        print("surface brightness output file will be created")
    if options['phot'] != None:
        params.flagphot=True
        print("output photometry file will be created")
    if options['minlevel'] != None:
        params.flagminlevel=True
    if options['sectors'] != None:
        params.flagsectors=True
    if options['object'] != None:
        params.flagobj=True
    if options['filter'] != None:
        params.flagband=True
    if options['snr'] != None:
        params.flagsnr=True
    if options['checkimg'] != None:
        params.flagcheck=True
    if options['noned'] != None:
        params.flagned=True
    if options['distmod'] != None:
        params.flagmod=True
    if options['magcor'] != None:
        params.flagmag=True
    if options['scalekpc'] != None:
        params.flagscale=True
    if options['sbdim'] != None:
        params.flagdim=True
    if options['model'] != None:
        params.flagmodel=True
        print("input model image will be used")
    if options['sky'] != None:
        params.flagsky=True
    if options['keep'] != None:
        params.flagkeep=True
    if options['ned'] != None:
        params.flagnedfile=True
    if options['gradsky'] != None:
        params.flagradsky=True
    if options['randsky'] != None:
        params.flagrandboxsky=True

    if options['skyRad'] != None:
        params.flagskyRad=True
    if options['skyRadmax'] != None:
        params.flagskyRadmax=True
    if options['skywidth'] != None:
        params.flagskywidth=True
    if options['skybox'] != None:
        params.flagskybox=True
    if options['skynum'] != None:
        params.flagskynum=True

    if options['distmax'] != None:
        params.flagdistmax=True

    if options['fwhm'] != None:
        params.flagfwhm=True




    # check for unrecognized options:
    sysopts=argv[2:]
    for idx,key in enumerate(sysopts):
        if not(key in OptionHandleList): 
            if key[0] == "-":
                print("WARNING: {} option not recognized ".format(key)) 

    if options['help'] != None:
        Help()

    ################## search arguments after the option:
    if params.flagpa == True:
        opt={}
        OptionHandle="-pa"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.parg=np.float(opt['pa'])

    if params.flagq == True:
        opt={}
        OptionHandle="-q"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.qarg=np.float(opt['q'])

    if params.flagranx[0] == True:
        opt={}
        OptionHandle="-ranx"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]

        params.rangex=opt["ranx"]
        if "-" in params.rangex:
            params.flagranx[1] = True
            params.ranx=opt['ranx']
        else:
            params.ranx=np.float(opt['ranx'])

    if params.flagrany[0]== True:
        opt={}
        OptionHandle="-rany"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]

        params.rangey=opt["rany"]
        if "-" in params.rangey:
            params.flagrany[1] = True
            params.rany=opt['rany']
        else:
            params.rany=np.float(opt['rany'])

    if params.flagdpi == True:
        opt={}
        OptionHandle="-dpi"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.dpival=np.int(opt['dpi'])

    if params.flagnoplot == True:
        params.dplot=False

    if params.flagminlevel == True:
        opt={}
        OptionHandle="-minlevel"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.minlevel=np.float(opt['minlevel'])

    if params.flagsectors == True:
        opt={}
        OptionHandle="-sectors"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.sectors=np.int(opt['sectors'])

    if params.flagobj == True:
        opt={}
        OptionHandle="-object"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.objname=np.str(opt['object'])

    if params.flagband == True:
        opt={}
        OptionHandle="-filter"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.band=np.str(opt['filter'])


    if params.flagmod == True:
        opt={}
        OptionHandle="-distmod"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.InDistMod=np.float(opt['distmod'])


    if params.flagmag == True:
        opt={}
        OptionHandle="-magcor"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.InMagCor=np.float(opt['magcor'])


    if params.flagscale == True:
        opt={}
        OptionHandle="-scalekpc"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.InScale=np.float(opt['scalekpc'])


    if params.flagdim == True:
        opt={}
        OptionHandle="-sbdim"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.InSbDim=np.float(opt['sbdim'])

    if params.flagmodel == True:
        opt={}
        OptionHandle="-model"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.inputmodel=np.str(opt['model'])


    if params.flagsky == True:
        opt={}
        OptionHandle="-sky"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.insky=np.float(opt['sky'])


    if params.flagnedfile == True:
        opt={}
        OptionHandle="-ned"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.nedfile=np.str(opt['ned'])


    if params.flagskyRad == True:
        opt={}
        OptionHandle="-skyRad"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.skyRad=np.float(opt['skyRad'])

    if params.flagskyRadmax == True:
        opt={}
        OptionHandle="-skyRadmax"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.skyRadmax=np.float(opt['skyRadmax'])


    if params.flagskybox == True:
        opt={}
        OptionHandle="-skybox"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.skybox=np.int(opt['skybox'])

    if params.flagskynum == True:
        opt={}
        OptionHandle="-skynum"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.skynum=np.int(opt['skynum'])


    if params.flagskywidth == True:
        opt={}
        OptionHandle="-skywidth"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.skywidth=np.int(opt['skywidth'])


    if params.flagdistmax == True:
        opt={}
        OptionHandle="-distmax"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.distmax=np.float(opt['distmax'])


    if params.flagfwhm == True:
        opt={}
        OptionHandle="-fwhm"
        opt[OptionHandle[1:]] = argv[argv.index(OptionHandle)+1]
        params.fwhm=np.float(opt['fwhm'])




    params.galfile= argv[1]

    return params



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



