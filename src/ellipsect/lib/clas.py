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



## class for parameters
class EllipSectConfig:

    #flags
    flaglogx=False
    flagq=False
    flagpa=False
    flagcomp=False
    flagpix=False
    flagranx=False
    flagrany=False
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
    #flagchi = False

    flagcheck=False
    flagned=False

    flagmod=False
    flagmag=False
    flagscale=False
    flagdim=False
   
    flagmodel=False

    flagsky=False
    flagkeep=True
    flagalax=False

    flagnedfile=False

    flagradsky=False
    flagrandboxsky=False


    flagskyRad = False
    flagskyRadmax = False
    flagskybox =  False
    flagskynum = False
    flagskywidth = False
    
    flagdistmax = False

    flagfwhm = False

    flagrep = False
    flagr90p = False
    flagr95p = False
    #flagcent = False 

    flagrmsky=True



    #init
    xc = 1 
    yc = 1 
    inxc=1    #same as above xc but used for input image
    inyc=1    #same as above yc but used for input image
 
    qarg=1
    parg=0
    ranx=1
    rany=1
    dplot=True

    dpival=100

    minlevel=0
    sectors=19

    title = False

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
   

    namefile = "none"
    namepng = "none.png"
    namesec = "none-gal.png"
    namemod = "none-mod.png"
    namemul = "none-mul.png"
    namesub = "none-sub.fits"
    namesig = "none-sig.fits"
    namesnr = "none-snr.fits"
    namechi = "none-chi.fits"
    namened = "none-ned.xml"
    namecheck = "none-check.fits"
    namering = "none-ring.fits"
    nameringmask = "none-ringmask.fits"

    namecube = "none-cub.png"


    nedfile = "default.xml"


    # sky parameters:
    skylevel = 0


    skyRad = 50 # minimum radius
    skybox = 20
    skynum = 20
    skywidth = 20

    fwhm = 2

    rep = 0
    r90p = 0
    r95p = 0

    distmax = 10

    brightness = 0 
    contrast = 1 

    cmap="viridis"


    numcomp = 1
    tot = 0 # total number of components


    # for computed gradsky
    gradskymean = 0
    gradskystd = 0
    gradskymed = 0

    randskymean = 0
    randskystd = 0
    randskymed = 0

    Aext = 0  # surface brightness correction for plots

    hconst = 67.8 
    omegam  = 0.308
    omegav = 0.692



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
    namechi="none-chi.fits"
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
    Aext = 0

    hconst = 67.8 
    omegam  = 0.308
    omegav = 0.692


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

    #from galfitcomps free parameters
    Active = np.array([])            #activate component  for galaxy

    # store the flags related to parameters
    PosXFree = np.array([])            #1)   
    PosYFree = np.array([])            #2)   
    MagFree = np.array([])             #3)
    RadFree = np.array([])             #4)
    ExpFree = np.array([])             #5)
    Exp2Free = np.array([])            #6)  for moffat
    Exp3Free = np.array([])            #7)  for moffat
                                   #8)  There is No 8 in any galfit model
    AxRatFree = np.array([])           #9)  AxisRatio
    PosAngFree = np.array([])          #10) position angle


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

    KronRad=np.array([])
    PetRad=np.array([])



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
    BICres = 0



class DataMge:
    ''' class for sector_photometry'''
    rad = 0 

    count = 0 
    angle = 0 
    anrad = 0 
    
    sb = 0 


class DataMgeComp:
    ''' class for sector_photometry per component'''
    rad = []

    count = []
    angle = []
    anrad = []
    
    sb = []

    sector = []


class DataTidal:
    ''' class for tidal function and other variables '''

    tidal = 0
    objchinu =0 
    bump = 0
    snr = 0
    stdsnr = 0 
    totsnr = 0
    rss = 0
    ndof = 0
    magalaper = 99  
    magmodaper = 99 


    totFlux = 0
    totMag = 99  
    BulgeToTotal = 0

    AbsMag = 99
    AbsMag2 = 99
    Lum = 0
    AICrit   = 0
    BICrit = 0
    BICres = 0

    EffRad = 0
    EffRad9 = 0
    EffRad3 = 0

    aell = 0
    bell = 0

class DataNed:
    ''' class for NED info '''

    GalExt = 0
    DistMod = 0
    DistMod2 = 0
    Scalekpc = 0 
    SbDim = 0



##### end of classes

### Dictionary for Absolute mag of the Sun taken from Willmer 2018
SunMag = {
        "U":5.61,"B":5.44,"V": 4.81,"R":4.43,"I":4.1,
        "J": 3.67,"H": 3.32,"K": 3.27,
        "u":5.49,"g":5.23,"r": 4.53,"i":4.19,"z":4.01,
        "L":3.26
        } 

####


