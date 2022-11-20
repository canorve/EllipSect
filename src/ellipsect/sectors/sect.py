
from ellipsect.lib.libs import *

from ellipsect import *

from ellipsect.lib.clas import GalfitParams
from ellipsect.lib.clas import GalfitComps 
from ellipsect.lib.clas import PhotAPI 

from ellipsect.inout.read import ReadGALFITout 
#from ellipsect.inout.read import GetWCS
from ellipsect.inout.read import ReadNComp 

from ellipsect.inout.plots import ShowCube 

from ellipsect.sectors.ellip import EllipSectors
from ellipsect.sectors.ellip import MulEllipSectors

from ellipsect.phot.phot import OutPhot

from ellipsect.lib.clas import EllipSectConfig


def SectorsGalfit(args):


    #note: remove this function
    ellconf = PassArgs(args) # from now on, ellconf is used instead of args

    #note: change this class for two class: one for the galfit header
    # and other for the galfit components


    #class for GALFIT's parameters
    galpar=GalfitParams()

    #class for GALFIT's components
    galcomps=GalfitComps()

    #note: change for one class 
    #class for output photometry 
    photapi=PhotAPI()

    #note: make one class for the ellipsect configuration

    ######################################
    ####### Read Galfit File #############
    #note: use two: one for the header and other for the components
    # create a class for reading
    ReadGALFITout(ellconf,galpar)
    ######################################
    ######################################

    if ellconf.flagq == True:
        galpar.q=ellconf.qarg

    if ellconf.flagpa == True:
        galpar.ang=ellconf.parg


    if ellconf.flagsky:
        galpar.skylevel=ellconf.insky


    #note: make one function to print all the configuration:
    str = "q = {} ".format(galpar.q)
    print(str)

    str = "pa = {} ".format(galpar.ang)
    print(str)

    ##
    str = "number of sectors = {}  ".format(ellconf.sectors)
    print(str)

    print("\nother parameters: \n")

    str = "Mag zeropoint = {} ".format(galpar.mgzpt)
    print(str)

    str = "Plate Scale = {} ".format(galpar.scale)
    print(str)

    str = "sky = {} ".format(galpar.skylevel)
    print(str)

    ##
    str = "minlevel = {} ".format(ellconf.minlevel)
    print(str)


    ##
    str = "for plots dpi = {} ".format(ellconf.dpival)
    print(str)
    ##

    print("grey angle at lower-bottom in multi-plot is measured from the galaxy's major axis ")
    print("red angle at upper-right in multi-plot is measured from the Y-axis (same as GALFIT)\n")


    print("In multi-plot, each color represents the same as the ones in the single-plot's legend")

    #names of the output files based on prefix of galfit output


    #note: make one function to save all the names: 
    root_ext = os.path.splitext(galpar.outimage)

    ellconf.namefile = root_ext[0]

    # names for the different png

    ellconf.namepng = ellconf.namefile + ".png"
    ellconf.namesec = ellconf.namefile + "-gal.png"
    ellconf.namemod = ellconf.namefile + "-mod.png"
    ellconf.namemul = ellconf.namefile + "-mul.png"
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



    if ellconf.flagsbout == True: 

        if not os.path.exists("sbfiles"):
            print("Creating directory for  output photometry ... ")
            os.makedirs("sbfiles")

        msg="prefix for surface brightness output file: {} ".format(ellconf.sboutput)
        print(msg)

    if ellconf.flagphot == True: 
        msg="output photometry file: {} ".format(ellconf.output)
        print(msg)

    #note: refactor this function move it above below reading the galfit header
    # read all the object components of the model in galfit.XX
    ReadNComp(ellconf.galfile,galpar.xc,galpar.yc,galcomps,ellconf.distmax)
    print("Number of components = ",len(galcomps.N))

    #note: move the ifs below to a function. Consider to create a class for reading
    #reading galaxy and model images from file
    errmsg="file {} does not exist".format(galpar.outimage)

    assert os.path.isfile(galpar.outimage), errmsg

    if ellconf.flagmodel == False:
        # hdu 1 => image   hdu 2 => model
        hdu = fits.open(galpar.outimage)
        galpar.img = (hdu[1].data.copy()).astype(float)
        galpar.model = (hdu[2].data.copy()).astype(float)
        galpar.imres = (hdu[3].data.copy()).astype(float)
        hdu.close()

    else:
        hdu = fits.open(galpar.inputimage)

        if galpar.flagidx:
            if galpar.flagnum:
                galpar.img = (hdu[galpar.imgidx,galpar.num].data).astype(float)
            else:    
                galpar.img = (hdu[galpar.imgidx].data).astype(float)
        else:
            galpar.img = (hdu[0].data).astype(float)

        hdu.close()

        hdu = fits.open(ellconf.inputmodel)
        galpar.model = (hdu[0].data).astype(float)
        hdu.close()

        galpar.imres = galpar.img - galpar.model

    # removing background from galaxy and model images 
    galpar.img = galpar.img - galpar.skylevel
    galpar.model = galpar.model - galpar.skylevel


    ### reading mask image from file

    if galpar.tempmask is not None:

        errmsg="file {} does not exist".format(galpar.tempmask)
        assert os.path.isfile(galpar.tempmask), errmsg

        if ellconf.flagmodel == False:
            hdu = fits.open(galpar.tempmask)
            mask = hdu[0].data
            galpar.mask=np.array(mask,dtype=bool)
            hdu.close()
        else:
            hdu = fits.open(galpar.maskimage)
            mask = hdu[0].data
            galpar.mask=np.array(mask,dtype=bool)
            hdu.close()

    else:
        galpar.mask=None
    ####
    #note: move to a function?
    ######################
    #shows the image cube#
    ######################{

    linewidth = 1.2

    if (ellconf.flagcomp):

        ell = Comp2Ellip(galpar, galcomps, linewidth)
    else:
        ell=[]


    #wcs = GetWCS(galpar.outimage) #removed

    ShowCube(galpar.outimage, namepng = ellconf.namecube, dpival 
             = ellconf.dpival, bri = ellconf.brightness, con = ellconf.contrast, 
             frac = ellconf.frac, fracmax = ellconf.fracmax,  cmap = ellconf.cmap, ellipse = ell)

    if ellconf.dplot:
        plt.pause(1.5)
 
    plt.close()


    #}
    #####################
    #####################



    #   numsectors=19
    #   numsectors=15
    numsectors=ellconf.sectors

    # minlevel=-100  # minimun value for sky
    # minlevel=15  # minimun value for sky
    minlevel=ellconf.minlevel  # minimun value for sky

    # initial values for image matrixes 
    sectgalax=sectmodel=sectcomps=[]

    #call to sectors_photometry for galaxy and model
    #note: divide the function below in two separated:
    sectgalax,sectmodel=SectPhot(galpar, ellconf, n_sectors=numsectors, minlevel=minlevel)

    
    if ellconf.flagcomp:
        #Note: sectors photometry for components always finished 
        # in minlevel = 0 regardless of the input -minlevel
        #sectcomps=SectPhotComp(galpar, ellconf, galcomps, n_sectors=numsectors, minlevel=minlevel)
        sectcomps=SectPhotComp(galpar, ellconf, galcomps, n_sectors=numsectors, minlevel=0)


    print("creating plots..")

    limx,limy=EllipSectors(ellconf, galpar, galcomps, sectgalax,sectmodel, sectcomps,n_sectors=numsectors)

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

    #note separate here in one for computation and other for plotting
    MulEllipSectors(ellconf, galpar, galcomps, sectgalax, sectmodel, sectcomps)


    print("multi-plot file: ", ellconf.namemul)


    if ellconf.dplot:
        plt.pause(1.5)

    plt.savefig(ellconf.namemul,dpi=ellconf.dpival)
    plt.close()


    ########################################################
    ############ Computing output photometry ###############
    ########################################################


    if ellconf.flagphot:
        print("Computing output photometry ... ")

        OutPhot(ellconf, galpar, galcomps, sectgalax, sectmodel, sectcomps, photapi)


    if galpar.tempmask != None:
        os.remove(galpar.tempmask) # removing temp mask file


    #note evalue how to eliminate the most class and variables
    PassVars(photapi,ellconf,galpar,galcomps)    

    return photapi





def SectPhot(galpar, ellconf, n_sectors=19, minlevel=0):
    """ calls to function sectors_photometry for galaxy and model """


    maskb=galpar.mask


    eps=1-galpar.q

    if ellconf.dplot:
        plt.clf()
        print("")

    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    if ellconf.flagmodel == False:
        yctemp=galpar.xc
        xctemp=galpar.yc
    else:
        yctemp=galpar.inxc
        xctemp=galpar.inyc


    # and angle is different as well:
    angsec=90-galpar.ang
    #    angsec=ang


    ###############################
    #  galaxy:


    sectgalax = sectors_photometry(galpar.img, eps, angsec, xctemp, yctemp, minlevel=minlevel,
            plot=ellconf.dplot, badpixels=maskb, n_sectors=n_sectors)


    if ellconf.dplot:
        plt.pause(1)  # Allow plot to appear on the screen
        plt.savefig(ellconf.namesec)

    ###################################################
    
    #  model: 
    # user input minlevel
    #sectmodel = sectors_photometry(galpar.model, eps, angsec, xctemp, yctemp,minlevel=minlevel,
    #        plot=ellconf.dplot, badpixels=maskb, n_sectors=n_sectors)
    # minlevel =0
    sectmodel = sectors_photometry(galpar.model, eps, angsec, xctemp, yctemp,minlevel=0,
            plot=ellconf.dplot, badpixels=maskb, n_sectors=n_sectors)


    if ellconf.dplot:
        plt.pause(1)  # Allow plot to appear on the screen
        plt.savefig(ellconf.namemod)



    return sectgalax,sectmodel


#sectors/sect.py
def SectPhotComp(galpar, ellconf, galcomps, n_sectors=19, minlevel=0):
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

    subimgs=[]

    mac=platform.system()

    if mac == 'Darwin':
        initcomp=1
    else:
        initcomp=2 #old galfit version
        initcomp=1 #new galfit version init subcomponents in 1

    cnt=0  # image =0 do not count
    while(cnt<len(galcomps.Comps)):
        if galcomps.Comps[cnt] == True:
            img = hdu[cnt+initcomp].data.astype(float)
            subimgs.append(img)
        cnt=cnt+1
    hdu.close()


    maskb=galpar.mask


    eps=1-galpar.q

    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    yctemp=galpar.xc
    xctemp=galpar.yc
    # and angle is different as well:
    angsec=90-galpar.ang

    epsmul=eps
    angsecmul=angsec
    #print("eps,angle mul ",epsmul,angsecmul)
    #############
    sectcomps=[]
    #sectmulcomps=[]
    n=0

    while(n<len(galcomps.N)):

        subim=subimgs[n]

        eps=1-galcomps.AxRat[n]
        angsec=90-galcomps.PosAng[n]


        if ellconf.flagcheck:
            scmp = sectors_photometry(subim, eps, angsec, xctemp, yctemp,minlevel=minlevel,plot=1, badpixels=maskb, n_sectors=n_sectors)
            plt.savefig("Comp"+str(n)+".png")
        else:
            scmp = sectors_photometry(subim, eps, angsec, xctemp, yctemp,minlevel=minlevel,plot=0, badpixels=maskb, n_sectors=n_sectors)


        #scmpmul = sectors_photometry(subim, epsmul, angsecmul, xctemp, yctemp,minlevel=minlevel,plot=0, badpixels=maskb, n_sectors=n_sectors)
        #plt.savefig("Cmul"+str(n)+".png")

        sectcomps.append(scmp)

        #sectmulcomps.append(scmpmul)

        n=n+1


    return sectcomps


def Comp2Ellip(galpar,galcomps,lw=1):
    ''' converts galfit component parameter into an Ellipse object''' 


    ellipses = [] 
    #col = 'red'

    N=len(galcomps.N)

    #color value
    values = range(N)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)



    
    for idx, item in enumerate(galcomps.N):


        # correcting coordinates
        xc=galcomps.PosX[idx] - galpar.xmin + 1
        yc=galcomps.PosY[idx] - galpar.ymin + 1


        pa = galcomps.PosAng[idx] + 90

        w = galcomps.Rad[idx]
        h = galcomps.Rad[idx]*galcomps.AxRat[idx]


        colorVal = scalarMap.to_rgba(values[idx])


        ell=Ellipse((xc, yc), width=w, height=h,angle=pa,
                     edgecolor=colorVal,
                     facecolor='none',
                     linewidth=lw)

        ellipses.append(ell)


    return ellipses




def PassVars(photapi,ellconf,galpar,galcomps):

    #################
    #from InputParams
    #ellconf=InputParams()

    #input file
    photapi.galfile=ellconf.galfile 

    #sb output file
    photapi.sboutput =ellconf.sboutput

    #output file
    photapi.output =ellconf.output

    # input image model
    photapi.inputmodel=ellconf.inputmodel


    # object name to search in NED
    photapi.objname=ellconf.objname
    photapi.namefile=ellconf.namefile
    photapi.namepng=ellconf.namepng
    photapi.namesec=ellconf.namesec
    photapi.namemod=ellconf.namemod
    photapi.namemul=ellconf.namemul
    photapi.namesub=ellconf.namesub
    photapi.namesig=ellconf.namesig
    photapi.namesnr=ellconf.namesnr
    photapi.namened=ellconf.namened
    photapi.namecheck=ellconf.namecheck
    photapi.namering=ellconf.namering
    photapi.nedfile=ellconf.nedfile


    #################
    #from GalfitParams
    #galpar=Galfitellconf()


    photapi.xc=galpar.xc
    photapi.yc=galpar.yc
    photapi.q=galpar.q
    photapi.ang=galpar.ang
    photapi.skylevel=galpar.skylevel
    photapi.scale=galpar.scale
    photapi.inputimage=galpar.inputimage
    photapi.outimage=galpar.outimage
    photapi.maskimage=galpar.maskimage
    photapi.mgzpt=galpar.mgzpt
    photapi.exptime=galpar.exptime
    photapi.tempmask=galpar.tempmask
    photapi.xmin=galpar.xmin
    photapi.xmax=galpar.xmax
    photapi.ymin=galpar.ymin
    photapi.ymax=galpar.ymax
    photapi.band=galpar.band
    photapi.inputimage=galpar.inputimage


    # from gradsky
    photapi.gradskymean =galpar.gradskymean
    photapi.gradskystd =galpar.gradskystd
    photapi.gradskymed =galpar.gradskymed

    # from randboxsky
    photapi.randskymean =galpar.randskymean
    photapi.randskystd =galpar.randskystd
    photapi.randskymed =galpar.randskymed

    #################
    # from GalfitComps
    #galcomps=GalfitComps()

    # init sub values
    # todos estos son arrays
    photapi.Comps=galcomps.Comps.copy()
    photapi.N=galcomps.N.copy()

    photapi.NameComp=galcomps.NameComp.copy()
    photapi.PosX=galcomps.PosX.copy()
    photapi.PosY=galcomps.PosY.copy()
    photapi.Mag=galcomps.Mag.copy()
    photapi.Rad=galcomps.Rad.copy()
    photapi.Exp=galcomps.Exp.copy()
    photapi.Exp2=galcomps.Exp2.copy()
    photapi.Exp3=galcomps.Exp3.copy()
                 
    photapi.AxRat=galcomps.AxRat.copy()
    photapi.PosAng =galcomps.PosAng.copy()
    photapi.skip=galcomps.skip.copy()
    photapi.freepar=galcomps.freepar.copy()

    # computed parameters:
    photapi.Rad50=galcomps.Rad50.copy()
    photapi.SerInd=galcomps.SerInd.copy()
    photapi.Rad50kpc=galcomps.Rad50kpc.copy()
    photapi.Rad50sec=galcomps.Rad50sec.copy()
    photapi.Rad90=galcomps.Rad90.copy()
    photapi.AbsMagComp=galcomps.AbsMag.copy()
    photapi.LumComp=galcomps.Lum.copy()
    photapi.Flux=galcomps.Flux.copy()
    photapi.PerLight=galcomps.PerLight.copy()
    photapi.me=galcomps.me.copy()
    photapi.mme=galcomps.mme.copy()
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
        ellconf.flagpa= True
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
        ellconf.flagsky= True
        ellconf.insky= args.sky

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



    return ellconf












