

from ellipsect.lib.libs import *



## class for parameters
class InputParams:

    #flags
    flaglogx=False
    flagq=False
    flagpa=False
    flagcomp=False
    flagpix=False
    flagranx=[False,False]
    flagrany=[False,False]
    flagnoplot=False
    flagrid=False
    flagdpi=False
    flagsbout=False
    flagphot=False
    flagminlevel=False
    flagsectors=False
    flagobj=False
    flagband=False
    flagweb=False # to check connection to ned
    flagsnr = False

    flagcheck=False
    flagned=False

    flagmod=False
    flagmag=False
    flagscale=False
    flagdim=False
   
    flagmodel=False

    flagsky=False
    flagkeep=False
    flagnedfile=False

    flagradsky=False
    flagrandboxsky=False


    flagskyRad = False
    flagskyRadmax = False
    flagskybox =  False
    flagskynum = False
    flagskywidth = False


    #init
    qarg=1
    parg=0
    ranx=1
    rany=1
    dplot=True

    dpival=100

    minlevel=0
    sectors=19



    #input file
    galfile= "galfit.01"

    #sb output file
    sboutput = "sbout.txt"

    #output file
    output = "out.txt"

    # input image model
    inputmodel="none.fits"

    # object name to search in NED

    objname="none"
    band="R"

    InDistMod=0
    InMagCor=0
    InScale=1
    InSbDim=0
   

    namefile="none"
    namepng="none.png"
    namesec="none-gal.png"
    namemod="none-mod.png"
    namemul="none-mul.png"
    namesub="none-sub.fits"
    namesig="none-sig.fits"
    namesnr="none-snr.fits"
    namened="none-ned.xml"
    namecheck="none-check.fits"
    namering="none-ring.fits"

    nedfile="default.xml"


    # sky parameters:
    insky=0

    skyRad = 50 # minimum radius
    skybox = 20
    skynum = 20
    skywidth = 20


#io/class.py
### class for Galfit parameters
class GalfitParams:

    xc=1        #for sectors_photometry
    yc=1        #for sectors_photometry
    inxc=1        #same as above xc but used for input image
    inyc=1        #same as above yc but used for input image
    q=1         #for sectors_photometry
    rad=0         #for sky rad 
    serind=0         #for sky rad 
    ang=0       #for sectors_photometry
    skylevel=0
    scale=1
    inputimage="galaxy.fits"
    outimage="galaxy-out.fits"
    maskimage="galaxy-mask.fits"
    mgzpt=25
    exptime=1
    tempmask="tempmask.fits"
    xmin=1
    xmax=2
    ymin=1
    ymax=2

    band="R"

    inputimage="galaxy.fits"

    imgidx="sci"
    flagidx=False
    num=1
    flagnum=False


    img = np.array([[1,1],[1,1]])
    model = np.array([[1,1],[1,1]])
    imres = np.array([[1,1],[1,1]])

    mask = np.array([[1,1],[1,1]])
    sigma = np.array([[1,1],[1,1]])
    imsnr = np.array([[1,1],[1,1]])


    # for computed gradsky
    gradskymean = 0
    gradskystd = 0
    gradskymed = 0

    randskymean = 0
    randskystd = 0
    randskymed = 0



#io/class.py
### class for Galfit components
class GalfitComps:

    # init sub values
    Comps=np.array([])  #remove this?
    N=np.array([])

    NameComp=np.array([])  #0)
    PosX=np.array([])            #1)   
    PosY=np.array([])            #2)   
    Mag=np.array([])             #3)
    Rad=np.array([])             #4)
    Exp=np.array([])             #5)
    Exp2=np.array([])            #6)  for moffat
    Exp3=np.array([])            #7)  for moffat
                                  #8)  There is No 8 in any galfit model
    AxRat=np.array([])           #9)  AxisRatio
    PosAng =np.array([])         #10) position angle
    skip=np.array([])            #z)  skip model
    freepar=np.array([])            # Number of free params

    # computed parameters:
    Rad50=np.array([])
    SerInd=np.array([])
    Rad50kpc=np.array([])
    Rad50sec=np.array([])
    Rad90=np.array([])
    AbsMag=np.array([])
    Lum=np.array([])
    Flux=np.array([])
    PerLight=np.array([])
    me=np.array([])
    mme=np.array([])
    kser = np.array([])


#class to comunicate externally:
class PhotAPI:

    #################
    #from InputParams

    #input file
    galfile= "galfit.01"

    #sb output file
    sboutput = "sbout.txt"

    #output file
    output = "out.txt"

    # input image model
    inputmodel="none.fits"


    # object name to search in NED
    objname="none"
    namefile="none"
    namepng="none.png"
    namesec="none-gal.png"
    namemod="none-mod.png"
    namemul="none-mul.png"
    namesub="none-sub.fits"
    namesig="none-sig.fits"
    namesnr="none-snr.fits"
    namened="none-ned.xml"
    namecheck="none-check.fits"
    namering="none-ring.fits"
    nedfile="default.xml"



    #################
    #from GalfitParams

    xc=1        #for sectors_photometry
    yc=1        #for sectors_photometry
    q=1         #for sectors_photometry
    ang=0       #for sectors_photometry
    skylevel=0
    scale=1
    inputimage="galaxy.fits"
    outimage="galaxy-out.fits"
    maskimage="galaxy-mask.fits"
    mgzpt=25
    exptime=1
    tempmask="tempmask.fits"
    xmin=1
    xmax=2
    ymin=1
    ymax=2
    band="R"
    inputimage="galaxy.fits"


    # from gradsky
    gradskymean = 0
    gradskystd = 0
    gradskymed = 0

    # from randboxsky
    randskymean = 0
    randskystd = 0
    randskymed = 0


    #################
    # from GalfitComps
    # init sub values
    Comps=np.array([])  #remove this?
    N=np.array([])

    NameComp=np.array([])  #0)
    PosX=np.array([])            #1)   
    PosY=np.array([])            #2)   
    Mag=np.array([])             #3)
    Rad=np.array([])             #4)
    Exp=np.array([])             #5)
    Exp2=np.array([])            #6)  for moffat
    Exp3=np.array([])            #7)  for moffat
                                  #8)  There is No 8 in any galfit model
    AxRat=np.array([])           #9)  AxisRatio
    PosAng =np.array([])         #10) position angle
    skip=np.array([])            #z)  skip model
    freepar=np.array([])            # Number of free params

    # computed parameters:
    Rad50=np.array([])
    SerInd=np.array([])
    Rad50kpc=np.array([])
    Rad50sec=np.array([])
    Rad90=np.array([])
    AbsMagComp=np.array([])
    LumComp=np.array([])
    Flux=np.array([])
    PerLight=np.array([])
    me=np.array([])
    mme=np.array([])
    kser = np.array([])




    #################
    # from output phot

    aell = 0
    bell = 0
    GalExt =0 
    DistMod=0
    DistMod2 = 0

    Scalekpc = 0
    SbDim = 0
    magalaper = 99
    magmodaper= 99 



    totFlux = 0
    totMag = 99  
    BulgeToTotal = 0
    tidal = 0
    objchinu = 0
    bump = 0
    snr = 0
    stdsnr = 0
    totsnr = 0
    rss   = 0
    ndof   = 0
    AbsMag = 99
    AbsMag2 = 99
    Lum = 0
    AICrit   = 0
    BICrit = 0


 
##### end of classes


