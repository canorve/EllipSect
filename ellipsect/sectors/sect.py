
from ellipsect.lib.libs import *

from ellipsect import *


from ellipsect.lib.clas import GalfitParams
from ellipsect.lib.clas import GalfitComps 
from ellipsect.lib.clas import PhotAPI 

from ellipsect.inout.read import ReadGALFITout 

from ellipsect.inout.read import ReadNComp 

from ellipsect.sectors.ellip import EllipSectors
from ellipsect.sectors.ellip import MulEllipSectors

from ellipsect.phot.phot import OutPhot



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
        initcomp=2

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

#/sectors/sect.py
def SectorsGalfit(params):


    print("angle in multi-plot is measured from the galaxy's major axis ")

    #class for GALFIT's parameters
    galpar=GalfitParams()

    #class for GALFIT's components
    galcomps=GalfitComps()

    #class for output photometry 
    photapi=PhotAPI()



    ######################################
    ####### Read Galfit File #############
    #  xc,yc,q,ang,skylevel,scale,file,mgzpt,exptime,mask=ReadGALFITout(params.galfile,galpars)
    ReadGALFITout(params.galfile,galpar)

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

    str = "sky = {} ".format(galpar.skylevel)
    print(str)

    ##
    str = "dpi = {} for plots ".format(params.dpival)
    print(str)
    ##

    ##
    str = "minlevel = {} ".format(params.minlevel)
    print(str)

    ##
    str = "number of sectors = {}  ".format(params.sectors)
    print(str)

    str = "Mag zeropoint = {} ".format(galpar.mgzpt)
    print(str)

    str = "Plate Scale = {} ".format(galpar.scale)
    print(str)



    (tmp)=galpar.outimage.split(".")

    params.namefile=tmp[0]

    # names for the different png

    params.namepng=params.namefile + ".png"
    params.namesec=params.namefile + "-gal.png"
    params.namemod=params.namefile + "-mod.png"
    params.namemul=params.namefile + "-mul.png"
    params.namesub=params.namefile + "-comp.fits"

    params.namesig=params.namefile + "-sig.fits"


    params.sboutput=params.namefile + "-sbout"
    params.output=params.namefile + "-out.txt"

    params.namened=params.namefile + "-ned.xml"



    params.namesnr=params.namefile + "-snr.fits"

    params.namecheck=params.namefile + "-check.fits"
    
    params.namering=params.namefile + "-ring.fits"



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


    if params.dplot:
        plt.pause(1.5)

    plt.savefig(params.namemul,dpi=params.dpival)
    plt.close()


    ########################################################
    ############ Computing output photometry ###############
    ########################################################

    PassVars(photapi,params,galpar,galcomps)    

    if params.flagphot:
        print("Computing output photometry ... ")

        OutPhot(params, galpar, galcomps, sectgalax, sectmodel, sectcomps, photapi)


    if galpar.tempmask != None:
        os.remove(galpar.tempmask) # removing temp mask file


    return photapi


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









