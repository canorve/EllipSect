
from ellipsect.lib.libs import *

from ellipsect import *

from ellipsect.lib.clas import GalfitParams
from ellipsect.lib.clas import GalfitComps 
from ellipsect.lib.clas import PhotAPI 

from ellipsect.inout.read import ReadGALFITout 
from ellipsect.inout.prt import printWelcome 

from ellipsect.inout.read import ReadNComp 
from ellipsect.inout.plots import ShowCube 

from ellipsect.sectors.ellip import EllipSectors
from ellipsect.sectors.ellip import MulEllipSectors

from ellipsect.phot.phot import OutPhot

from ellipsect.lib.clas import InputParams


def SectorsGalfit(args):


    params = PassArgs(args) # from now on, params is used instead of args


    printWelcome() # print version and description 

    #class for GALFIT's parameters
    galpar=GalfitParams()

    #class for GALFIT's components
    galcomps=GalfitComps()

    #class for output photometry 
    photapi=PhotAPI()


    ######################################
    ####### Read Galfit File #############
    ReadGALFITout(params,galpar)
    ######################################
    ######################################

    if params.flagq == True:
        galpar.q=params.qarg

    if params.flagpa == True:
        galpar.ang=params.parg


    if params.flagsky:
        galpar.skylevel=params.insky

    str = "q = {} ".format(galpar.q)
    print(str)

    str = "pa = {} ".format(galpar.ang)
    print(str)

    ##
    str = "number of sectors = {}  ".format(params.sectors)
    print(str)

    print("\nother parameters: \n")

    str = "Mag zeropoint = {} ".format(galpar.mgzpt)
    print(str)

    str = "Plate Scale = {} ".format(galpar.scale)
    print(str)

    str = "sky = {} ".format(galpar.skylevel)
    print(str)

    ##
    str = "minlevel = {} ".format(params.minlevel)
    print(str)


    ##
    str = "for plots dpi = {} ".format(params.dpival)
    print(str)
    ##

    print("grey angle at lower-bottom in multi-plot is measured from the galaxy's major axis ")
    print("red angle at upper-right in multi-plot is measured from the Y-axis (same as GALFIT)\n")


    print("In multi-plot, each color represents the same as the ones in the single-plot's legend")

    #names of the output files based on prefix of galfit output

    #(tmp)=galpar.outimage.split(".")
    root_ext = os.path.splitext(galpar.outimage)

    params.namefile = root_ext[0]

    # names for the different png

    params.namepng = params.namefile + ".png"
    params.namesec = params.namefile + "-gal.png"
    params.namemod = params.namefile + "-mod.png"
    params.namemul = params.namefile + "-mul.png"
    params.namesub = params.namefile + "-comp.fits"

    params.namesig = params.namefile + "-sig.fits"


    params.sboutput = params.namefile + "-sbout"
    params.output = params.namefile + "-out.txt"

    params.namened = params.namefile + "-ned.xml"



    params.namesnr = params.namefile + "-snr.fits"

    params.namecheck = params.namefile + "-check.fits"
    
    params.namering = params.namefile + "-ring.fits"
    
    params.nameringmask = params.namefile + "-ringmask.fits"

    params.namecube = params.namefile + "-cub.png"



    if params.flagsbout == True: 

        if not os.path.exists("sbfiles"):
            print("Creating directory for  output photometry ... ")
            os.makedirs("sbfiles")

        msg="prefix for surface brightness output file: {} ".format(params.sboutput)
        print(msg)

    if params.flagphot == True: 
        msg="output photometry file: {} ".format(params.output)
        print(msg)


    # read all the object components of the model in galfit.XX
    ReadNComp(params.galfile,galpar.xc,galpar.yc,galcomps,params.distmax)
    print("Number of components = ",len(galcomps.N))


    #reading galaxy and model images from file
    errmsg="file {} does not exist".format(galpar.outimage)

    assert os.path.isfile(galpar.outimage), errmsg

    if params.flagmodel == False:
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

        hdu = fits.open(params.inputmodel)
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

        if params.flagmodel == False:
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

    ######################
    #shows the image cube#
    ######################{

    linewidth = 1.2

    if (params.flagcomp):

        ell = Comp2Ellip(galpar,galcomps,linewidth)
    else:
        ell=[]


    ShowCube(galpar.outimage,namepng=params.namecube,dpival=params.dpival,frac=params.frac,cmap=params.cmap,ellipse=ell)

    if params.dplot:
        plt.pause(1.5)
 
    plt.close()


    #}
    #####################
    #####################



    #   numsectors=19
    #   numsectors=15
    numsectors=params.sectors

    # minlevel=-100  # minimun value for sky
    # minlevel=15  # minimun value for sky
    minlevel=params.minlevel  # minimun value for sky

    # initial values for image matrixes 
    sectgalax=sectmodel=sectcomps=[]

    #call to sectors_photometry for galaxy and model
    sectgalax,sectmodel=SectPhot(galpar, params, n_sectors=numsectors, minlevel=minlevel)

    
    if params.flagcomp:
        #Note: sectors photometry for components always finished 
        # in minlevel = 0 regardless of the input -minlevel
        #sectcomps=SectPhotComp(galpar, params, galcomps, n_sectors=numsectors, minlevel=minlevel)
        sectcomps=SectPhotComp(galpar, params, galcomps, n_sectors=numsectors, minlevel=0)



    print("creating plots..")

    limx,limy=EllipSectors(params, galpar, galcomps, sectgalax,sectmodel, sectcomps,n_sectors=numsectors)

    print("plot file: ", params.namepng)
 


    ##############################################
    ##############################################
    ##############################################

    if params.dplot:
        plt.pause(1.5)
    plt.savefig(params.namepng,dpi=params.dpival)
    #plt.close()

    ########################################################
    ################ Multiplots: ###########################
    ########################################################

    print("creating multi-plots..")

    MulEllipSectors(params, galpar, galcomps, sectgalax, sectmodel, sectcomps)


    print("multi-plot file: ", params.namemul)


    if params.dplot:
        plt.pause(1.5)

    plt.savefig(params.namemul,dpi=params.dpival)
    plt.close()


    ########################################################
    ############ Computing output photometry ###############
    ########################################################


    if params.flagphot:
        print("Computing output photometry ... ")

        OutPhot(params, galpar, galcomps, sectgalax, sectmodel, sectcomps, photapi)


    if galpar.tempmask != None:
        os.remove(galpar.tempmask) # removing temp mask file


    PassVars(photapi,params,galpar,galcomps)    

    return photapi





def SectPhot(galpar, params, n_sectors=19, minlevel=0):
    """ calls to function sectors_photometry for galaxy and model """


    maskb=galpar.mask


    eps=1-galpar.q

    if params.dplot:
        plt.clf()
        print("")

    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    if params.flagmodel == False:
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
            plot=params.dplot, badpixels=maskb, n_sectors=n_sectors)


    if params.dplot:
        plt.pause(1)  # Allow plot to appear on the screen
        plt.savefig(params.namesec)

    ###################################################
    
    #  model: 
    # user input minlevel
    #sectmodel = sectors_photometry(galpar.model, eps, angsec, xctemp, yctemp,minlevel=minlevel,
    #        plot=params.dplot, badpixels=maskb, n_sectors=n_sectors)
    # minlevel =0
    sectmodel = sectors_photometry(galpar.model, eps, angsec, xctemp, yctemp,minlevel=0,
            plot=params.dplot, badpixels=maskb, n_sectors=n_sectors)


    if params.dplot:
        plt.pause(1)  # Allow plot to appear on the screen
        plt.savefig(params.namemod)



    return sectgalax,sectmodel


#sectors/sect.py
def SectPhotComp(galpar, params, galcomps, n_sectors=19, minlevel=0):
    """ calls to function sectors_photometry for subcomponents """

    if (params.flagphot) and (not(os.path.isfile(params.namesig))):

        if ((os.path.isfile(params.namesub)) and (params.flagkeep)):
            print("using existing subcomponent model image file *-comp.fits")

            print("running galfit to create sigma image ...")

            rungal = "galfit -outsig {}".format(params.galfile)
            errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            # changing name to sigma image
            runchg = "mv sigma.fits {}".format(params.namesig)
            errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            errmsg="file {} does not exist".format(params.namesig)
            assert os.path.isfile(params.namesig), errmsg

        else:

            print("running galfit to create sigma image and individual model images...")

            rungal = "galfit -o3 -outsig {}".format(params.galfile)
            errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            # changing name to subcomponents
            runchg = "mv subcomps.fits {}".format(params.namesub)
            errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            errmsg="file {} does not exist".format(params.namesub)
            assert os.path.isfile(params.namesub), errmsg

            # changing name to sigma image

            runchg = "mv sigma.fits {}".format(params.namesig)
            errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            errmsg="file {} does not exist".format(params.namesig)
            assert os.path.isfile(params.namesig), errmsg

    else: 

        if ((os.path.isfile(params.namesub)) and (params.flagkeep)):
            print("using existing subcomponent model image file *-comp.fits")
        else:

            print("running galfit to create individual model images...")

            rungal = "galfit -o3 {}".format(params.galfile)
            errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            # changing name to subcomponents
            runchg = "mv subcomps.fits {}".format(params.namesub)
            errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

            errmsg="file {} does not exist".format(params.namesub)
            assert os.path.isfile(params.namesub), errmsg


        if (os.path.isfile(params.namesig) and (params.flagphot)):
            print("using existing sigma image")

    ##

    hdu = fits.open(params.namesub)

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


        if params.flagcheck:
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




def PassVars(photapi,params,galpar,galcomps):

    #################
    #from InputParams
    #params=InputParams()

    #input file
    photapi.galfile=params.galfile 

    #sb output file
    photapi.sboutput =params.sboutput

    #output file
    photapi.output =params.output

    # input image model
    photapi.inputmodel=params.inputmodel


    # object name to search in NED
    photapi.objname=params.objname
    photapi.namefile=params.namefile
    photapi.namepng=params.namepng
    photapi.namesec=params.namesec
    photapi.namemod=params.namemod
    photapi.namemul=params.namemul
    photapi.namesub=params.namesub
    photapi.namesig=params.namesig
    photapi.namesnr=params.namesnr
    photapi.namened=params.namened
    photapi.namecheck=params.namecheck
    photapi.namering=params.namering
    photapi.nedfile=params.nedfile


    #################
    #from GalfitParams
    #galpar=GalfitParams()


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
    '''function to pass arguments from args to params'''

    # Note: This function shouldn't exist, but since 
    # I didn't know about the argparse library when 
    # I start this project, I have to create this function
    # to pass the arguments from argparse to my old parsing arguments
    # this is the minimum thing to do without modifying the rest 
    # of the code.

    #class for user's parameters
    params=InputParams()

    # passing to params
    ##########################
    params.galfile= args.GalFile 
    ##########################

    #options without arguments

    if args.logx:
        params.flaglogx=True

    if args.comp:
        params.flagcomp=True
 
    if args.pix:
        params.flagpix=True
 
    if args.grid:
        params.flagrid=True
 
    if args.sbout:
        params.flagsbout=True
 
    if args.noplot:
        params.flagnoplot=True
        params.dplot=False

    if args.phot:
        params.flagphot=True
 
    if args.checkimg:
        params.flagcheck=True
 
    if args.noned:
        params.flagned=True
 
    if args.gradsky:
        params.flagradsky=True

    if args.randsky:
        params.flagrandboxsky=True

    if args.snr:
        params.flagsnr=True
 
    if args.keep:
        params.flagkeep=True

    if args.galax:
        params.flagalax=True

    if args.allskypx:
        params.flagrmsky=False



    #options with arguments

    if args.center:
        params.flagcent = True
        params.xc = args.center[0]
        params.yc = args.center[1]


    if args.axisrat:
        params.flagq = True
        params.qarg = args.axisrat  

    if args.posangle:
        params.flagpa= True
        params.parg = args.posangle  


    if args.ranx:
        params.flagranx=True
        params.ranx=args.ranx

    if args.rany:
        params.flagrany=True
        params.rany=args.rany


    if args.dotsinch:
        params.flagdpi = True
        params.dpival = args.dotsinch

    if args.minlevel:
        params.flagminlevel= True
        params.minlevel= args.minlevel

    if args.sectors:
        params.flagsectors= True
        params.sectors= args.sectors

    if args.object:
        params.flagobj= True
        params.objname= args.object

    if args.filter:
        params.flagband= True
        params.band= args.filter

    if args.distmod:
        params.flagmod = True
        params.InDistMod = args.distmod

    if args.magcor:
        params.flagmag= True
        params.InMagCor= args.magcor

    if args.scalekpc:
        params.flagscale= True
        params.InScale= args.scalekpc

    if args.sbdim:
        params.flagdim= True
        params.InSbDim= args.sbdim

    if args.model:
        params.flagmodel= True
        params.inputmodel= args.model

    if args.sky:
        params.flagsky= True
        params.insky= args.sky

    if args.ned:
        params.flagnedfile= True
        params.nedfile= args.ned

    if args.radinit:
        params.flagskyRad= True
        params.skyRad= args.radinit

    if args.skyradmax:
        params.flagskyRadmax= True
        params.skyRadmax= args.skyradmax

    if args.skybox:
        params.flagskybox= True
        params.skybox= args.skybox
  
    if args.skynum:
        params.flagskynum= True
        params.skynum= args.skynum

    if args.skywidth:
        params.flagskywidth= True
        params.skywidth= args.skywidth

    if args.distmax:
        params.flagdistmax= True
        params.distmax= args.distmax
   
    if args.fwhm:
        params.flagfwhm= True
        params.fwhm= args.fwhm

    if args.frac:
        params.frac = args.frac

    if args.cmap:
        params.cmap = args.cmap



    return params












