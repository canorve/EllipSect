

from ellipsect.lib.libs import *


from ellipsect.sectors.num import GetK
from ellipsect.phot.tidal import Tidal
from ellipsect.phot.ned import NED


from ellipsect.sectors.num import Re90

#phot/phot.py
### Dictionary for Absolute mag of the Sun taken from Willmer 2018
SunMag = {
        "U":5.61,"B":5.44,"V": 4.81,"R":4.43,"I":4.1,
        "J": 3.67,"H": 3.32,"K": 3.27,
        "u":5.49,"g":5.23,"r": 4.53,"i":4.19,"z":4.01,
        "L":3.26
        } 

####





def OutPhot(params, galpar, galcomps, sectgalax, sectmodel, sectcomps, photapi):
    """ Output photometry for further analysis """


    # masks to identify components: 

    maskmag=(galcomps.NameComp != "ferrer") & (galcomps.NameComp != "nuker") & (galcomps.NameComp != "edgedisk") & (galcomps.NameComp != "king") 
    maskgalax = (galcomps.NameComp == "sersic") | (galcomps.NameComp == "devauc") | (galcomps.NameComp == "expdisk")  | (galcomps.NameComp == "gaussian") 


    maskdisk = (galcomps.NameComp == "expdisk") # ignore this at the moment: #| (galcomps.NameComp == "edgedisk") 
    maskbulge = (galcomps.NameComp == "sersic") | (galcomps.NameComp == "devauc") | (galcomps.NameComp == "moffat") | (galcomps.NameComp == "ferrer") | (galcomps.NameComp == "king") | (galcomps.NameComp == "gaussian") | (galcomps.NameComp == "psf")
    masksersic = (galcomps.NameComp == "sersic") 
    maskexp= (galcomps.NameComp == "expdisk")
    maskdevauc=(galcomps.NameComp == "devauc")
    maskgauss=(galcomps.NameComp == "gaussian")


    if maskmag.any():

        galcomps.Flux[maskmag]=10**((galpar.mgzpt - galcomps.Mag[maskmag])/2.5)
        totFlux=galcomps.Flux[maskmag].sum()
        totMag=-2.5*np.log10(totFlux) + galpar.mgzpt
        #print("total magnitud = ",totMag)
        #print("total Flux = ",totFlux)
    #else:
    #    print("Total magnitud can not be computed with the actual galfit functions")



    if maskdisk.any():
        BulgeFlux = 10**((galpar.mgzpt -  galcomps.Mag[maskbulge])/2.5)
        DiskFlux  = 10**((galpar.mgzpt -  galcomps.Mag[maskdisk])/2.5)
        totBulgeF = BulgeFlux.sum()
        totDiskF =  DiskFlux.sum()
        BulgeToTotal= totBulgeF / (totBulgeF + totDiskF)
        #print("BulgeToTotal = ",BulgeToTotal)
    else:
        BulgeToTotal = 1
        #print("BulgeToTotal = 1")



    ################## Computing component variables ###########

    # sersic index
    galcomps.SerInd[masksersic] = galcomps.Exp[masksersic] 
    galcomps.SerInd[maskdevauc] = 4
    galcomps.SerInd[maskexp]    = 1
    galcomps.SerInd[maskgauss]  = 0.5


    # effective radius
    galcomps.Rad50[masksersic] = galcomps.Rad[masksersic] 
    galcomps.Rad50[maskdevauc] = galcomps.Rad[maskdevauc] 

    galcomps.Rad50[maskexp] = 1.678*galcomps.Rad[maskexp] 
    galcomps.Rad50[maskgauss] = 0.5*galcomps.Rad[maskgauss] 

    ################
    # computing Rad 90% light with aproximation taken from my Thesis

    galcomps.Rad50sec[maskgalax] = galcomps.Rad50[maskgalax] * galpar.scale 

    #galcomps.Rad90[maskgalax] = galcomps.Rad50[maskgalax] * (1.53 + 0.73 * galcomps.SerInd[maskgalax] + 0.07 * galcomps.SerInd[maskgalax]**2) 
    galcomps.Rad90[maskgalax] = Re90(galcomps.Rad50[maskgalax],galcomps.SerInd[maskgalax] )

    ######

    print("Rad90 is the radius at 90% of total light  ")



    galcomps.kser[maskgalax]  = GetK(galcomps.SerInd[maskgalax])

    # computing meanme and me 
    galcomps.mme[maskgalax]  = galcomps.Mag[maskgalax]  + 2.5 * np.log10(2 * np.pi * galcomps.AxRat[maskgalax] * galcomps.Rad50sec[maskgalax]**2 )


    fn = (( galcomps.AxRat[maskgalax] * galcomps.SerInd[maskgalax] * np.exp( galcomps.kser[maskgalax])) / (galcomps.kser[maskgalax] ** (2 * galcomps.SerInd[maskgalax] )) ) * ( np.exp(scipy.special.gammaln(2*galcomps.SerInd[maskgalax])) )

    galcomps.me[maskgalax] = galcomps.mme[maskgalax] +  2.5 * np.log10( fn )



    Num=len(galcomps.Flux[maskmag])


    if (Num > 0):

        galcomps.PerLight[maskmag]= (galcomps.Flux[maskmag] / totFlux )
        namecomp=galcomps.NameComp[maskmag]
        N=galcomps.N[maskmag]

        N=N.astype(int)

    #n=0
    #while(n<Num):

    #    print("Num, componente, %perlight ",N[n],namecomp[n],galcomps.PerLight[n])

    #    n+=1

   
    #print("Number of free params: ",int(galcomps.freepar.sum()))


    # create sigma image
    if(not(os.path.isfile(params.namesig))):

        print("running galfit to create sigma image...")

        rungal = "galfit -o2 -outsig {} ".format(params.galfile) 
        errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        runchg = "mv sigma.fits {}".format(params.namesig)
        errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        errmsg="file {} does not exist".format(params.namesig)
        assert os.path.isfile(params.namesig), errmsg
    else:
        if not(params.flagcomp):
            print("using existing sigma image ")


    # ignore warnings from Card too long
    if not sys.warnoptions:
        warnings.simplefilter("ignore")


    # call to Tidal
    (tidal,objchinu,bump,snr,stdsnr,totsnr,rss,ndof,magalaper,magmodaper)=Tidal(params, galpar, galcomps, sectgalax, 2)

    # returns warnings to normal
    if not sys.warnoptions:
        warnings.simplefilter("default")




    print("galaxy mag using sectors_photometry aperture = {:.3f} ".format(magalaper))
    print("Model mag using sectors_photometry aperture = {:.3f}".format(magmodaper))


    #print("Tidal = ",tidal)  
    #print("Local Chinu = ",objchinu)  
    #print("Bumpiness = ",bump)  
    #print("mean SNR = ",snr)  
    #print("std SNR = ",stdsnr)  
    #print("total SNR sum over area = ",totsnr)  
    #print("Residual sum of squares =  ",rss)  
    #print("degrees of freedom = ",ndof)  
   
    header = fits.getheader(galpar.inputimage)
    #hdu = fits.open(galpar.inputimage)
    #header = hdu[0].header
    #hdu.close()




    if not(params.flagobj):
        if "OBJECT" in header: 
            params.objname=header["OBJECT"] 
            params.flagobj=True
            print("using object name: {} to search in NED ".format(params.objname))            

        elif "TARGNAME" in header:  
            params.objname=header["TARGNAME"] 
            params.flagobj=True
            print("using object name: {} to search in NED ".format(params.objname))            
        else:
            print("WARNING: name for object not found in header nor it was provided by user") 
            print("Luminosity and absolute magnitude will not be computed") 
            print("Luminosity and absolute magnitude for individual components will have wrong quantities ") 

    else:
        print("using object name: {} to search in NED ".format(params.objname))            


    if not(params.flagband):
        if "BAND" in header: 
            params.flagband = True
            params.band=header["BAND"] 
        elif "FILTER" in header: 
            params.flagband = True
            params.band=header["FILTER"] 
        elif "FILTNAM" in header: 
            params.flagband = True
            params.band=header["FILTNAM"] 
        elif "FILTNAM1" in header: 
            params.flagband = True
            params.band=header["FILTNAM1"] 
        else:
            print("WARNING: filter not found. using default filter: ",params.band) 
            print("use -filter option to change band") 
    else:
        print("using {} band to correct for galactic extinction ".format(params.band)) 

    if (params.flagobj and not(params.flagned)): 
        (GalExt,DistMod,DistMod2,Scalekpc,SbDim)=NED(params, galpar, galcomps)
    else:
        if params.flagned: 
            print("No search in NED. Lum and abs Mag will not be computed")

        GalExt=0
        DistMod=0
        DistMod2=0
        Scalekpc=0
        SbDim=0


    if maskmag.any():


        CorMag = totMag - GalExt # corrected by galactic extinction 
            
        AbsMag=CorMag - DistMod # No K correction applied

        AbsMag2=CorMag - DistMod2 # No K correction applied

        # per component: 
        CompCorMag = galcomps.Mag[maskmag] - GalExt # corrected by galactic extinction 
        #galcomps.AbsMag[maskmag] = CompCorMag - DistMod # No K correction applied
        #No K correction applied:
        galcomps.AbsMag[maskmag] = CompCorMag - DistMod2 # AbsMag using distance modulus independent of z 

        if params.band in SunMag:
            MSun = SunMag[params.band]

            #Lum = 10**((MSun - AbsMag)/2.5) Luminosity will  now be computed using AbsMag2
            Lum = 10**((MSun - AbsMag2)/2.5)
            # per component 
            galcomps.Lum[maskmag]= 10**((MSun - galcomps.AbsMag[maskmag])/2.5)

        else:
            print("Absolute Magnitude for Band {} was not found. Check filter name ".format(params.band))
            print("Luminosity will not be computed.")

            Lum = 0

        #print("Magnitud Absoluta",AbsMag)
        #print("Magnitud Absoluta using Distance Modulus independen of z ",AbsMag2)
        #print("check references in ",params.namened)


    if maskgalax.any():

        #if (params.flagweb and params.flagobj and not(params.flagned)):

        galcomps.Rad50kpc[maskgalax] = galcomps.Rad50[maskgalax] * galpar.scale * Scalekpc

        galcomps.mme[maskgalax] = galcomps.mme[maskgalax] - GalExt - SbDim
        galcomps.me[maskgalax] = galcomps.me[maskgalax] - GalExt - SbDim
 
        #else:
        #    print("mean surface brightness at Re is not corrected for galactic extintion nor surface brightness dimming ")





    ################  INFORMATION CRITERIA ####################


    freepar=int(galcomps.freepar.sum())

    npix = ndof + freepar

    #    ;  AKAIKE INFORMATION CRITERION
    #    ; AIC = chi^2 + 2k

    AICrit = objchinu * ndof + 2*freepar
    

    #    ; BAYESIAN INFORMATION CRITERION
    #    ; BIC = chi^2 + k * ln(n)

    BICrit = objchinu * ndof + freepar * np.log(npix)

    ## for output only: 
    stidxg = np.argsort(sectgalax.radius)

    mgerad=sectgalax.radius[stidxg]

    aell = mgerad.max() 
    bell = mgerad.max() * galpar.q



    #######  file output:  ######

    print("Creating output photometry file: ",params.output)

    OUTPHOT = open (params.output,"w")

    lineout= "#   Output photometry for {} file \n".format(galpar.outimage)
    OUTPHOT.write(lineout)

    lineout= "#\n"
    OUTPHOT.write(lineout)

    lineout = "# Photometry for object: {} (use -object option to change it) \n".format(params.objname)
    OUTPHOT.write(lineout)

    lineout = "# All photometric quantities are computed for filter {} (use -filter option to change it) \n".format(params.band)
    OUTPHOT.write(lineout)

    lineout = "# Magnitudes are not corrected by K-Correction \n"
    OUTPHOT.write(lineout)

    lineout= "# Total magnitude and other variables do NOT include: \n"
    OUTPHOT.write(lineout)

    lineout= "# ferrer, nuker, edgedisk and king components.  \n"
    OUTPHOT.write(lineout)

    lineout = "# Some photometric variables are computed within an ellipse defined by sectors_photometry \n"
    OUTPHOT.write(lineout)

    lineout = "# This ellipse has axis a = {:.2f} and b = {:.2f} centered at xc, yc \n".format(aell,bell)
    OUTPHOT.write(lineout)

    lineout= "#\n"
    OUTPHOT.write(lineout)

    lineout= "# sectors_photometry was used with q={}, pa={} and minlevel = {} \n".format(galpar.q,galpar.ang,params.minlevel)
    OUTPHOT.write(lineout)

    lineout= "# OutImage = {}  MgZpt = {}  \n".format(galpar.outimage,galpar.mgzpt)
    OUTPHOT.write(lineout)

    lineout= "# exptime = {}  plate scale = {} ''/pix \n".format(galpar.exptime,galpar.scale)
    OUTPHOT.write(lineout)

    lineout= "# xc = {:.2f}  yc = {:.2f}  sky = {:.2f} \n".format(galpar.xc, galpar.yc, galpar.skylevel)
    OUTPHOT.write(lineout)
    
    lineout = "# Gal. Extinction = {}; Distance Mod. = {}; Distance Mod. (z independent) = {} \n".format(GalExt,DistMod,DistMod2)
    OUTPHOT.write(lineout)

    lineout = "# cosmology corrected scale = {} kpc/arcsec; Surface brightness dimming (mag) = {} \n".format(Scalekpc,SbDim)
    OUTPHOT.write(lineout)

    lineout = "# Check references in NED file {} \n\n".format(params.namened)
    OUTPHOT.write(lineout)


    lineout = "total apparent mag of the galaxy (using ellipse aperture on image) = {:.3f}  \n".format(magalaper)
    OUTPHOT.write(lineout)

    lineout = "total apparent mag of the model (using ellipse aperture on image) = {:.3f}  \n".format(magmodaper)
    OUTPHOT.write(lineout)


    if maskmag.any():
        #Flux=10**((galpar.mgzpt -  galcomps.Mag[maskmag])/2.5)
        #totFlux=Flux.sum()
        #totMag=-2.5*np.log10(totFlux) + galpar.mgzpt

        lineout = "total apparent mag from model parameters (without corrections) = {:.3f}  \n".format(totMag)
        OUTPHOT.write(lineout)

        lineout = "total flux (without corrections) = {:.3f}  \n".format(totFlux)
        OUTPHOT.write(lineout)
    else:
        lineout="Total magnitude can not be computed with the actual galfit functions \n"
        OUTPHOT.write(lineout)


    lineout="Bulge To Total Ratio = {:.3f} \n".format(BulgeToTotal)
    OUTPHOT.write(lineout)

    lineout="Tidal = {:.3f} \n".format(tidal)
    OUTPHOT.write(lineout)

    lineout="Local Chinu = {:.3f} \n".format(objchinu)
    OUTPHOT.write(lineout)

    lineout="Bumpiness = {:.3f} \n".format(bump)
    OUTPHOT.write(lineout)

    lineout = "mean SNR = {:.3f} \n".format(snr)
    OUTPHOT.write(lineout)

    lineout = "std SNR = {:.3f} \n".format(stdsnr)
    OUTPHOT.write(lineout)

    lineout = "total SNR sum = {:.3f} \n".format(totsnr)  
    OUTPHOT.write(lineout)

    lineout = "Residual sum of squares (RSS) = {:.3f}  \n".format(rss)  
    OUTPHOT.write(lineout)

    lineout= "degrees of freedom = {} \n".format(ndof)  
    OUTPHOT.write(lineout)


    lineout="Number of free params = {} \n".format(int(galcomps.freepar.sum()))
    OUTPHOT.write(lineout)
   

    if maskmag.any():


        if (params.flagweb and params.flagobj and not(params.flagned)):

            #CorMag = totMag - GalExt # corrected by galactic extinction 
            
            #AbsMag=CorMag - DistMod # No K correction applied

            #AbsMag2=CorMag - DistMod2 # No K correction applied

            lineout="Absolute Mag = {:.3f} \n".format(AbsMag)
            OUTPHOT.write(lineout)

            lineout="Absolute Mag using Dist Mod independent of z = {:.3f} \n".format(AbsMag2)
            OUTPHOT.write(lineout)

            lineout="Luminosity = {:.3f} (10^10 solar lum) using Dist Mod independent of z  \n".format(Lum/1e10)
            OUTPHOT.write(lineout)



    lineout= "Akaike Information Criterion = {:.3f} \n".format(AICrit)  
    OUTPHOT.write(lineout)

    lineout = "Bayesian Information Criterion = {:.3f} \n".format(BICrit)  
    OUTPHOT.write(lineout)

    lineout= "\n"
    OUTPHOT.write(lineout)

    ############################

    if params.flagradsky:

        mean = galpar.gradskymean
        std = galpar.gradskystd
        median = galpar.gradskymed


        lineout="computed sky with the gradient method:  mean = {:.2f} , std = {:.2f}, median = {} \n".format(mean,std,median)

        OUTPHOT.write(lineout)
 
    if params.flagrandboxsky:

        mean = galpar.randskymean
        std = galpar.randskystd
        median = galpar.randskymed

        lineout="computed sky with the random box method:  mean = {:.2f} , std = {:.2f}, median = {} \n".format(mean,std,median)
        OUTPHOT.write(lineout)


    lineout = "\n"  
    OUTPHOT.write(lineout)

    #############################

    lineout = "#########################################\n"  
    OUTPHOT.write(lineout)

    lineout = "# Photometric properties per component: #\n"  
    OUTPHOT.write(lineout)

    lineout = "#########################################\n"  
    OUTPHOT.write(lineout)

    lineout= "\n"
    OUTPHOT.write(lineout)


    lineout = "# Number   Component   FractionLight    me        <me>        AbsMag      Luminosity      Rad90       Re   \n"  
    OUTPHOT.write(lineout)

    lineout = "#                                    (mag/'')   (mag/'')      (mag)    (10^10 SolarLum)   (pix)      (kpc)   \n"  
    OUTPHOT.write(lineout)


    for idx, item in enumerate(galcomps.N) :
        lineout= "    {0:^2} {1:^17} {2:^10.3f} {3:^10.3f} {4:^10.3f} {5:^14.3f} {6:^14.3f} {7:^10.3f} {8:^10.3f}  \n".format(galcomps.N[idx],galcomps.NameComp[idx],galcomps.PerLight[idx],galcomps.me[idx],galcomps.mme[idx],galcomps.AbsMag[idx],galcomps.Lum[idx]/1e10,galcomps.Rad90[idx],galcomps.Rad50kpc[idx])
        OUTPHOT.write(lineout)

    OUTPHOT.close()

    # save variables for output class

    photapi.aell =aell
    photapi.bell =bell 
    photapi.GalExt =GalExt 
    photapi.DistMod=DistMod
    photapi.DistMod2 = DistMod2

    photapi.Scalekpc =Scalekpc
    photapi.SbDim = SbDim 
    photapi.magalaper =magalaper
    photapi.magmodaper=magmodaper 

    photapi.totFlux = totFlux 
    photapi.totMag = totMag 
    photapi.BulgeToTotal =BulgeToTotal
    photapi.tidal = tidal 
    photapi.objchinu = objchinu
    photapi.bump = bump 
    photapi.snr = snr
    photapi.stdsnr =stdsnr 
    photapi.totsnr =totsnr 
    photapi.rss   = rss 
    photapi.ndof   = ndof 
    photapi.AbsMag = AbsMag 
    photapi.AbsMag2 = AbsMag2 
    photapi.Lum = Lum 
    photapi.AICrit   = AICrit 
    photapi.BICrit = BICrit


