
from ellipsect.lib.libs import *

from ellipsect import *


from ellipsect.inout.read import ReadGALFITout 
from ellipsect.inout.read import prefixNames

from ellipsect.inout.plots import plotCube 

from ellipsect.sectors.ellip import EllipSectors
from ellipsect.sectors.ellip import MulEllipSectors

from ellipsect.phot.phot import OutPhot

from ellipsect.lib.clas import EllipSectConfig


from ellipsect.inout.galfit  import Galfit 
from ellipsect.inout.galfit  import numComps
from ellipsect.inout.galfit  import readDataImg

from ellipsect.inout.prt import printEllinfo

from ellipsect.sky.sky import skyCall


from ellipsect.lib.clas import PhotAPI 


def SectorsGalfit(args):


    ellconf = PassArgs(args) # from now on, ellconf is used instead of args


    ######################################
    ####### Read Galfit File #############


    galhead = Galfit.ReadHead(ellconf.galfile)
    galcomps = Galfit.ReadComps(ellconf.galfile)
    galsky = Galfit.ReadSky(ellconf.galfile)


    if ellconf.flagsky == False:

        ellconf.skylevel = galsky.sky


    ReadGALFITout(ellconf, galhead, galcomps) 
    ######################################
    ######################################

    printEllinfo(ellconf, galhead) #print parameter info

    N = numComps(galcomps)
    print("Number of components = ",N)

    ellconf.tot = N

    #creates names of the output files based on prefix of galfit output
    prefixNames(ellconf, galhead.outimage)



    if ellconf.flagsbout == True: 

        if not os.path.exists("sbfiles"):
            print("Creating directory for output photometry ... ")
            os.makedirs("sbfiles")

        msg="prefix for surface brightness output file: {} ".format(ellconf.sboutput)
        print(msg)

    if ellconf.flagphot == True: 
        msg="output photometry file: {} ".format(ellconf.output)
        print(msg)



    dataimg = readDataImg(ellconf, galhead)

    # removing background from galaxy and model images 
    dataimg.img = dataimg.img - ellconf.skylevel
    dataimg.model = dataimg.model - ellconf.skylevel


    plotCube(ellconf, galhead, galcomps) #plots the cube image


    #   numsectors=19
    #   numsectors=15
    numsectors = ellconf.sectors

    # minlevel=-100  # minimun value for sky
    # minlevel=15  # minimun value for sky
    minlevel = ellconf.minlevel  # minimun value for sky

    #call to sectors_photometry for galaxy and model

    sectgalax = SectPhot(ellconf, dataimg, n_sectors = numsectors, minlevel = minlevel, fit='gal')
    sectmodel = SectPhot(ellconf, dataimg, n_sectors = numsectors, minlevel = minlevel, fit='mod' )

    
    sectcomps=[]
    if ellconf.flagcomp:
        #Note: sectors photometry for components always finished 
        # in minlevel = 0 regardless of the input -minlevel

        sectcomps = SectPhotComp(ellconf, dataimg, galcomps, n_sectors = numsectors, minlevel = 0)
        

    #computing sky
    skyCall(ellconf, galhead, galcomps)



    print("creating plots..")

    limx, limy = EllipSectors(ellconf, galhead, galcomps, sectgalax, sectmodel, sectcomps, n_sectors = numsectors)

    print("plot file: ", ellconf.namepng)
 


    ##############################################
    ##############################################
    ##############################################

    if ellconf.dplot:
        plt.pause(1.5)
    plt.savefig(ellconf.namepng,dpi=ellconf.dpival)
    #plt.close()

    ########################################################
    ################ Multiplots: ###########################
    ########################################################



    print("creating multi-plots..")

    MulEllipSectors(ellconf, galhead, galcomps, sectgalax, sectmodel, sectcomps)


    print("multi-plot file: ", ellconf.namemul)


    if ellconf.dplot:
        plt.pause(1.5)

    plt.savefig(ellconf.namemul,dpi=ellconf.dpival)
    plt.close()


    ########################################################
    ############ Computing output photometry ###############
    ########################################################


    # save variables for output class
    photapi = PhotAPI() 


    if ellconf.flagphot:
        print("Computing output photometry ... ")

        photapi = OutPhot(ellconf, dataimg, galhead, galcomps, sectgalax, sectmodel, sectcomps, photapi)
        


    if galhead.tempmask != None:
        os.remove(galhead.tempmask) # removing temp mask file


    PassVars(photapi, ellconf, galhead, galcomps)    


    return photapi





def SectPhot(ellconf, dataimg, n_sectors = 19, minlevel = 0, fit = 'gal'):
    """ calls to function sectors_photometry for galaxy and model """


    maskb = dataimg.mask

    eps = 1 - ellconf.qarg

    if ellconf.dplot:
        plt.clf()
        print("")

    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    if ellconf.flagmodel == False:
        yctemp = ellconf.xc
        xctemp = ellconf.yc
    else:
        yctemp = ellconf.inxc
        xctemp = ellconf.inyc


    # and angle is different as well:
    angsec = 90 - ellconf.parg
    #    angsec=ang


    ###################################################
    if fit == 'gal':
    #  galaxy:

        sectimg = sectors_photometry(dataimg.img, eps, angsec, xctemp, yctemp, minlevel = minlevel,
                plot = ellconf.dplot, badpixels = maskb, n_sectors = n_sectors)


        if ellconf.dplot:
            plt.pause(1)  # Allow plot to appear on the screen
            plt.savefig(ellconf.namesec)


    if fit == 'mod':
    #  model: 
        sectimg = sectors_photometry(dataimg.model, eps, angsec, xctemp, yctemp,minlevel=0,
                plot = ellconf.dplot, badpixels = maskb, n_sectors = n_sectors)


        if ellconf.dplot:
            plt.pause(1)  # Allow plot to appear on the screen
            plt.savefig(ellconf.namemod)



    return sectimg


#sectors/sect.py
def SectPhotComp(ellconf, dataimg, galcomps, n_sectors=19, minlevel=0):
    """ calls to function sectors_photometry for subcomponents """

    if (ellconf.flagphot) and (not(os.path.isfile(ellconf.namesig))):

        if ((os.path.isfile(ellconf.namesub)) and (ellconf.flagkeep)):
            print("using existing subcomponent model image file *-comp.fits")

            print("running galfit to create sigma image ...")

            rungal = "galfit -outsig {}".format(ellconf.galfile)
            errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            # changing name to sigma image
            runchg = "mv sigma.fits {}".format(ellconf.namesig)
            errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            errmsg="file {} does not exist".format(ellconf.namesig)
            assert os.path.isfile(ellconf.namesig), errmsg

        else:

            print("running galfit to create sigma image and individual model images...")

            rungal = "galfit -o3 -outsig {}".format(ellconf.galfile)
            errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            # changing name to subcomponents
            runchg = "mv subcomps.fits {}".format(ellconf.namesub)
            errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            errmsg="file {} does not exist".format(ellconf.namesub)
            assert os.path.isfile(ellconf.namesub), errmsg

            # changing name to sigma image

            runchg = "mv sigma.fits {}".format(ellconf.namesig)
            errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            errmsg="file {} does not exist".format(ellconf.namesig)
            assert os.path.isfile(ellconf.namesig), errmsg

    else: 

        if ((os.path.isfile(ellconf.namesub)) and (ellconf.flagkeep)):
            print("using existing subcomponent model image file *-comp.fits")
        else:

            print("running galfit to create individual model images...")

            rungal = "galfit -o3 {}".format(ellconf.galfile)
            errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            # changing name to subcomponents
            runchg = "mv subcomps.fits {}".format(ellconf.namesub)
            errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            errmsg="file {} does not exist".format(ellconf.namesub)
            assert os.path.isfile(ellconf.namesub), errmsg


        if (os.path.isfile(ellconf.namesig) and (ellconf.flagphot)):
            print("using existing sigma image")

    ##

    hdu = fits.open(ellconf.namesub)

    subimgs = []

    mac = platform.system()

    if mac == 'Darwin':  #probably remove this part of the code
        initcomp = 1
    else:
        initcomp = 2 #old galfit version
        initcomp = 1 #new galfit version init subcomponents in 1

    cnt=0  # image = 0 do not count


    while(cnt < len(galcomps.N)):
        if galcomps.Active[cnt] == True:
            img = hdu[cnt + initcomp].data.astype(float)
            subimgs.append(img)
        cnt = cnt + 1
    hdu.close()


    maskb = dataimg.mask


    eps = 1 - ellconf.qarg

    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    yctemp = ellconf.xc
    xctemp = ellconf.yc
    # and angle is different as well:
    angsec = 90 - ellconf.parg

    epsmul = eps
    angsecmul = angsec
    #############
    sectcomps=[]
    n = 0


    #masksel = galcomps.Active == True 

    i = 0

    while(n < len(galcomps.N)):


        subim = subimgs[i]

        if galcomps.Active[n] == True:


            eps = 1 - galcomps.AxRat[n]
            angsec = 90 - galcomps.PosAng[n]


            if ellconf.flagcheck:

                scmp = sectors_photometry(subim, eps, angsec, xctemp, yctemp, minlevel = minlevel, 
                                            plot = 1, badpixels = maskb, n_sectors = n_sectors)

                plt.savefig("Comp"+str(n)+".png")

            else:
                scmp = sectors_photometry(subim, eps, angsec, xctemp, yctemp, minlevel = minlevel, 
                                            plot = 0, badpixels = maskb, n_sectors = n_sectors)



            sectcomps.append(scmp)

            i = i + 1

        n = n + 1

    return sectcomps



def PassVars(photapi, ellconf, galhead, galcomps):

    #################
    #from InputParams
    #ellconf=InputParams()

    #input file
    photapi.galfile = ellconf.galfile 

    #sb output file
    photapi.sboutput = ellconf.sboutput

    #output file
    photapi.output = ellconf.output

    # input image model
    photapi.inputmodel = ellconf.inputmodel


    # object name to search in NED
    photapi.objname = ellconf.objname
    photapi.namefile = ellconf.namefile
    photapi.namepng = ellconf.namepng
    photapi.namesec = ellconf.namesec
    photapi.namemod = ellconf.namemod
    photapi.namemul = ellconf.namemul
    photapi.namesub = ellconf.namesub
    photapi.namesig = ellconf.namesig
    photapi.namesnr = ellconf.namesnr
    photapi.namened = ellconf.namened
    photapi.namecheck = ellconf.namecheck
    photapi.namering = ellconf.namering
    photapi.nedfile = ellconf.nedfile
    photapi.band = ellconf.band

    #################
    #from Galhead and ellconf


    photapi.xc = ellconf.xc
    photapi.yc = ellconf.yc
    photapi.q = ellconf.qarg
    photapi.ang = ellconf.parg
    photapi.skylevel = ellconf.skylevel
    photapi.scale = galhead.scale
    photapi.inputimage = galhead.inputimage
    photapi.outimage = galhead.outimage
    photapi.maskimage = galhead.maskimage
    photapi.mgzpt = galhead.mgzpt
    photapi.exptime = galhead.exptime
    photapi.tempmask = galhead.tempmask
    photapi.xmin = galhead.xmin
    photapi.xmax = galhead.xmax
    photapi.ymin = galhead.ymin
    photapi.ymax = galhead.ymax



    # from gradsky
    photapi.gradskymean = ellconf.gradskymean
    photapi.gradskystd = ellconf.gradskystd
    photapi.gradskymed = ellconf.gradskymed

    # from randboxsky
    photapi.randskymean = ellconf.randskymean
    photapi.randskystd = ellconf.randskystd
    photapi.randskymed = ellconf.randskymed

    #################
    # from GalfitComps
    #galcomps=GalfitComps()

    # init sub values
    #photapi.Comps=galcomps.Comps.copy()
    photapi.N = galcomps.N.copy()

    photapi.NameComp = galcomps.NameComp.copy()
    photapi.PosX = galcomps.PosX.copy()
    photapi.PosY = galcomps.PosY.copy()
    photapi.Mag = galcomps.Mag.copy()
    photapi.Rad = galcomps.Rad.copy()
    photapi.Exp = galcomps.Exp.copy()
    photapi.Exp2 = galcomps.Exp2.copy()
    photapi.Exp3 = galcomps.Exp3.copy()
                 
    photapi.AxRat = galcomps.AxRat.copy()
    photapi.PosAng = galcomps.PosAng.copy()
    photapi.skip = galcomps.skip.copy()

    photapi.Active = galcomps.Active.copy()    

    # store the flags related to parameters
    photapi.PosXFree = galcomps.PosXFree.copy()
    photapi.PosYFree = galcomps.PosYFree.copy()
    photapi.MagFree = galcomps.MagFree.copy()
    photapi.RadFree = galcomps.RadFree.copy()
    photapi.ExpFree = galcomps.ExpFree.copy()
    photapi.Exp2Free = galcomps.Exp2Free.copy()
    photapi.Exp3Free = galcomps.Exp3Free.copy()
    photapi.AxRatFree = galcomps.AxRatFree.copy()
    photapi.PosAngFree = galcomps.PosAngFree.copy()




    # computed parameters:
    photapi.Rad50 = galcomps.Rad50.copy()
    photapi.SerInd = galcomps.SerInd.copy()
    photapi.Rad50kpc = galcomps.Rad50kpc.copy()
    photapi.Rad50sec = galcomps.Rad50sec.copy()
    photapi.Rad90 = galcomps.Rad90.copy()
    photapi.AbsMagComp = galcomps.AbsMag.copy()
    photapi.LumComp = galcomps.Lum.copy()
    photapi.Flux = galcomps.Flux.copy()
    photapi.PerLight = galcomps.PerLight.copy()
    photapi.me = galcomps.me.copy()
    photapi.mme = galcomps.mme.copy()
    photapi.kser = galcomps.kser.copy()




def PassArgs(args):
    '''function to pass arguments from args to ellconf'''

    # Note: This function shouldn't exist, but since 
    # I didn't know about the argparse library when 
    # I started this project, I have to create this function
    # to pass the arguments from argparse to my old parsing arguments
    # this is the minimum thing to do without modifying the rest 
    # of the code.

    #class for user's parameters
    ellconf = EllipSectConfig()

    # passing to ellconf
    ##########################
    ellconf.galfile = args.GalFile 
    ##########################

    #options without arguments

    if args.logx:
        ellconf.flaglogx=True

    if args.comp:
        ellconf.flagcomp=True
 
    if args.pix:
        ellconf.flagpix=True
 
    if args.grid:
        ellconf.flagrid=True
 
    if args.sbout:
        ellconf.flagsbout=True
 
    if args.noplot:
        ellconf.flagnoplot=True
        ellconf.dplot=False

    if args.phot:
        ellconf.flagphot=True
 
    if args.checkimg:
        ellconf.flagcheck=True
 
    if args.noned:
        ellconf.flagned=True
 
    if args.gradsky:
        ellconf.flagradsky=True

    if args.randsky:
        ellconf.flagrandboxsky=True

    if args.snr:
        ellconf.flagsnr=True
 
    if args.keep:
        ellconf.flagkeep=True

    if args.galax:
        ellconf.flagalax=True

    if args.allskypx:
        ellconf.flagrmsky=False



    #options with arguments

    if args.center:
        ellconf.flagcent = True
        ellconf.xc = args.center[0]
        ellconf.yc = args.center[1]


    if args.axisrat:
        ellconf.flagq = True
        ellconf.qarg = args.axisrat  

    if args.posangle:
        ellconf.flagpa = True
        ellconf.parg = args.posangle  


    if args.ranx:
        ellconf.flagranx=True
        ellconf.ranx=args.ranx

    if args.rany:
        ellconf.flagrany=True
        ellconf.rany=args.rany


    if args.dotsinch:
        ellconf.flagdpi = True
        ellconf.dpival = args.dotsinch

    if args.minlevel:
        ellconf.flagminlevel= True
        ellconf.minlevel= args.minlevel

    if args.sectors:
        ellconf.flagsectors= True
        ellconf.sectors= args.sectors

    if args.object:
        ellconf.flagobj= True
        ellconf.objname= args.object

    if args.filter:
        ellconf.flagband= True
        ellconf.band= args.filter

    if args.distmod:
        ellconf.flagmod = True
        ellconf.InDistMod = args.distmod

    if args.magcor:
        ellconf.flagmag= True
        ellconf.InMagCor= args.magcor

    if args.scalekpc:
        ellconf.flagscale= True
        ellconf.InScale= args.scalekpc

    if args.sbdim:
        ellconf.flagdim= True
        ellconf.InSbDim= args.sbdim

    if args.model:
        ellconf.flagmodel= True
        ellconf.inputmodel= args.model

    if args.sky:
        ellconf.flagsky = True
        ellconf.skylevel = args.sky

    if args.ned:
        ellconf.flagnedfile= True
        ellconf.nedfile= args.ned

    if args.radinit:
        ellconf.flagskyRad= True
        ellconf.skyRad= args.radinit

    if args.skyradmax:
        ellconf.flagskyRadmax= True
        ellconf.skyRadmax= args.skyradmax

    if args.skybox:
        ellconf.flagskybox= True
        ellconf.skybox= args.skybox
  
    if args.skynum:
        ellconf.flagskynum= True
        ellconf.skynum= args.skynum

    if args.skywidth:
        ellconf.flagskywidth= True
        ellconf.skywidth= args.skywidth

    if args.distmax:
        ellconf.flagdistmax= True
        ellconf.distmax= args.distmax
   
    if args.fwhm:
        ellconf.flagfwhm= True
        ellconf.fwhm= args.fwhm

    if args.brightness:
        ellconf.brightness = args.brightness


    if args.contrast:
        ellconf.contrast = args.contrast


    if args.frac:
        ellconf.frac = args.frac

    if args.fracmax:
        ellconf.fracmax = args.fracmax




    if args.cmap:
        ellconf.cmap = args.cmap

    if args.numcomp:
        ellconf.numcomp = args.numcomp




    return ellconf












