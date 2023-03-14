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



from ellipsect.lib.libs import  *


from ellipsect.sectors.num import GetK
from ellipsect.phot.tidal import Tidal
from ellipsect.phot.ned import NED


from ellipsect.sectors.num import Re90
from ellipsect.sectors.num import RadGamma

from ellipsect.sectors.num import KronRadius
from ellipsect.sectors.num import solvePet

from ellipsect.lib.clas import DataTidal
from ellipsect.lib.clas import SunMag
from ellipsect.lib.clas import DataNed

from ellipsect.inout.prt import printPhot


from ellipsect.inout.galfit  import numParFree
from ellipsect.inout.galfit  import GetReff 



#phot/phot.py



def OutPhot(ellconf, dataimg, galhead, galcomps, sectgalax, sectmodel, sectcomps,photapi):
    """ Computes Output photometry """





    # masks to identify components: 
    maskmag = ((galcomps.NameComp != "ferrer") & (galcomps.NameComp != "nuker") 
                & (galcomps.NameComp != "edgedisk") & (galcomps.NameComp != "king") 
                & (galcomps.Active == True))

    maskgalax = (((galcomps.NameComp == "sersic") | (galcomps.NameComp == "devauc") | 
                (galcomps.NameComp == "expdisk")  | (galcomps.NameComp == "gaussian")) & 
                (galcomps.Active == True))


    maskdisk = (galcomps.NameComp == "expdisk") & (galcomps.Active == True)
    # ignore this at the moment: #| (galcomps.NameComp == "edgedisk") 
    maskbulge = (((galcomps.NameComp == "sersic") | (galcomps.NameComp == "devauc") | 
                (galcomps.NameComp == "moffat") | (galcomps.NameComp == "ferrer") | 
                (galcomps.NameComp == "king") | (galcomps.NameComp == "gaussian") | 
                (galcomps.NameComp == "psf")) & (galcomps.Active == True))

    masksersic = (galcomps.NameComp == "sersic")  & (galcomps.Active == True)
    maskexp = (galcomps.NameComp == "expdisk") & (galcomps.Active == True)
    maskdevauc = (galcomps.NameComp == "devauc") & (galcomps.Active == True)
    maskgauss = (galcomps.NameComp == "gaussian") & (galcomps.Active == True)
    masknuker = (galcomps.NameComp == "nuker") & (galcomps.Active == True)

    datatidal = DataTidal()


    if maskmag.any():

        galcomps.Flux[maskmag] = 10**((galhead.mgzpt - galcomps.Mag[maskmag])/2.5)
        datatidal.totFlux = galcomps.Flux[maskmag].sum()
        datatidal.totMag = -2.5*np.log10(datatidal.totFlux) + galhead.mgzpt

    else:

        datatidal.totFlux = 0
        datatidal.totMag = 99


    if maskdisk.any():

        BulgeFlux = 10**((galhead.mgzpt -  galcomps.Mag[maskbulge])/2.5)
        DiskFlux  = 10**((galhead.mgzpt -  galcomps.Mag[maskdisk])/2.5)
        totBulgeF = BulgeFlux.sum()
        totDiskF =  DiskFlux.sum()
        datatidal.BulgeToTotal = totBulgeF / (totBulgeF + totDiskF)
        #print("BulgeToTotal = ",BulgeToTotal)
    else:
        datatidal.BulgeToTotal = 1
        #print("BulgeToTotal = 1")



    ################## Computing component variables ###########

    # sersic index
    galcomps.SerInd[masksersic] = galcomps.Exp[masksersic] 

    galcomps.SerInd[maskdevauc] = 4
    galcomps.Exp[maskdevauc] = 4
    
    galcomps.SerInd[maskexp]    = 1
    galcomps.Exp[maskexp]    = 1

    galcomps.SerInd[maskgauss]  = 0.5
    galcomps.Exp[maskgauss]  = 0.5


    # effective radius
    #check the name rad50 y Rad add it to the class?
    galcomps.Rad50[masksersic] = galcomps.Rad[masksersic] 
    galcomps.Rad50[maskdevauc] = galcomps.Rad[maskdevauc] 

    galcomps.Rad50[maskexp] = 1.678*galcomps.Rad[maskexp] 
    galcomps.Rad50[maskgauss] = 0.5*galcomps.Rad[maskgauss] 

    ################

    galcomps.Rad50sec[maskgalax] = galcomps.Rad50[maskgalax] * galhead.scale 

    #galcomps.Rad90[maskgalax] = galcomps.Rad50[maskgalax] * (1.53 + 0.73 * galcomps.SerInd[maskgalax] + 0.07 * galcomps.SerInd[maskgalax]**2) 
    galcomps.Rad90[maskgalax] = Re90(galcomps.Rad50[maskgalax], galcomps.SerInd[maskgalax] )
    print("Rad90 is the radius at 90% of total light of the component")

    galcomps.KronRad[maskgalax] = KronRadius(galcomps.Rad50[maskgalax], 
                                            galcomps.Rad50[maskgalax], galcomps.SerInd[maskgalax])


    galcomps.PetRad[maskgalax]= solvePet(galcomps.SerInd[maskgalax])
    print("Petrosian Radius is computed with 1/nu = 0.2 (in Re units)")

 
    if masknuker.any():
 
       
        # if the component is Nuker, galcomps.Rad90 represents  the gamma radius
        galcomps.Rad90[masknuker] = RadGamma(galcomps.Rad[masknuker], galcomps.Exp[masknuker], 
                                            galcomps.Exp2[masknuker],galcomps.Exp3[masknuker])


        print("For Nuker component, Rad90 is the gamma radius ")
    ######




    galcomps.kser[maskgalax]  = GetK(galcomps.SerInd[maskgalax])

    # computing meanme and me 
    galcomps.mme[maskgalax] = galcomps.Mag[maskgalax] + 2.5*np.log10(2*np.pi*galcomps.AxRat[maskgalax]* 
                                                                            galcomps.Rad50sec[maskgalax]**2)


    fn = (((galcomps.AxRat[maskgalax]*galcomps.SerInd[maskgalax]* 
            np.exp(galcomps.kser[maskgalax]))/(galcomps.kser[maskgalax]**
                                            (2*galcomps.SerInd[maskgalax])))*
            (np.exp(sc.gammaln(2*galcomps.SerInd[maskgalax]))))

    galcomps.me[maskgalax] = galcomps.mme[maskgalax] + 2.5*np.log10(fn)



    Num = len(galcomps.Flux[maskmag])


    if (Num > 0):

        galcomps.PerLight[maskmag] = (galcomps.Flux[maskmag]/datatidal.totFlux)
        namecomp = galcomps.NameComp[maskmag]
        N = galcomps.N[maskmag]

        N=N.astype(int)


    # create sigma image
    if(not(os.path.isfile(ellconf.namesig))):

        print("running galfit to create sigma image...")

        rungal = "galfit -o2 -outsig {} ".format(ellconf.galfile) 
        errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        runchg = "mv sigma.fits {}".format(ellconf.namesig)
        errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        errmsg="file {} does not exist".format(ellconf.namesig)
        assert os.path.isfile(ellconf.namesig), errmsg
    else:
        if not(ellconf.flagcomp):
            print("using existing sigma image ")


    # ignore warnings from Card too long
    if not sys.warnoptions:
        warnings.simplefilter("ignore")


    ################################################################
    ################################################################
    print('computing the effective radius for the whole galaxy')
    eff = 0.5
    EffRad, tempmag = GetReff().GetReSer(galhead, galcomps, eff, ellconf.parg)

    eff = 0.9
    EffRad9, tempmag = GetReff().GetReSer(galhead, galcomps, eff, ellconf.parg)

    eff = 0.3
    EffRad3, tempmag = GetReff().GetReSer(galhead, galcomps, eff, ellconf.parg)

    datatidal.EffRad = EffRad
    datatidal.EffRad9 = EffRad9
    datatidal.EffRad3 = EffRad3
    ################################################################
    ################################################################

    if maskmag.any():
        datatidal.aell = datatidal.EffRad9
        datatidal.bell = datatidal.EffRad9*ellconf.qarg

    else:

        print("model contains non-sersic functions.") 
        print("90% Radius from sector_photometry will be used instead")
        frac = .90
        datatidal.aell = mgerad.max()*frac 
        datatidal.bell = mgerad.max()*ellconf.qarg*frac



    # call to Tidal
    datatidal = Tidal(datatidal, ellconf, dataimg, galhead, galcomps, sectgalax, 2)


    # returns warnings msgs to normal
    if not sys.warnoptions:
        warnings.simplefilter("default")



    print("galaxy mag using sectors_photometry aperture = {:.3f}".format(datatidal.magalaper))
    print("Model mag using sectors_photometry aperture = {:.3f}".format(datatidal.magmodaper))

  
    header = fits.getheader(galhead.inputimage)


    if not(ellconf.flagobj):
        if "OBJECT" in header: 
            ellconf.objname=header["OBJECT"] 
            ellconf.flagobj=True
            print("using object name: {} to search in NED ".format(ellconf.objname))            

        elif "TARGNAME" in header:  
            ellconf.objname=header["TARGNAME"] 
            ellconf.flagobj=True
            print("using object name: {} to search in NED ".format(ellconf.objname))            
        else:
            print("WARNING: name for object not found in header nor it was provided by user") 
            print("Luminosity and absolute magnitude will not be computed") 
            print("Luminosity and absolute magnitude for individual components will have wrong quantities ") 

    else:
        print("using object name: {} to search in NED ".format(ellconf.objname))            


    if not(ellconf.flagband):
        if "BAND" in header: 
            ellconf.flagband = True
            ellconf.band=header["BAND"] 
        elif "FILTER" in header: 
            ellconf.flagband = True
            ellconf.band=header["FILTER"] 
        elif "FILTNAM" in header: 
            ellconf.flagband = True
            ellconf.band=header["FILTNAM"] 
        elif "FILTNAM1" in header: 
            ellconf.flagband = True
            ellconf.band=header["FILTNAM1"] 
        else:
            print("WARNING: filter not found. Using default filter: ",ellconf.band) 
            print("use --filter option to change band") 
    else:
        print("using {} band to correct for galactic extinction ".format(ellconf.band)) 


    if (ellconf.flagobj and not(ellconf.flagned)): 
        dataned = NED(ellconf)
    else:
        if ellconf.flagned: 
            print("No search in NED because it is indicated by the user")
            print("Lum and abs Mag will not be computed")
        elif not(ellconf.flagobj):
            print("No search in NED. Object name not provided or not found in Header")
            print("Lum and abs Mag will not be computed")

        dataned = DataNed()
        dataned.GalExt = 0
        dataned.DistMod = 0
        dataned.DistMod2 = 0
        dataned.Scalekpc = 0
        dataned.SbDim = 0


    if maskmag.any():


        CorMag = datatidal.totMag - dataned.GalExt # corrected by galactic extinction 
            
        datatidal.AbsMag=CorMag - dataned.DistMod # No K correction applied

        datatidal.AbsMag2=CorMag - dataned.DistMod2 # No K correction applied

        # per component: 
        CompCorMag = galcomps.Mag[maskmag] - dataned.GalExt # corrected by galactic extinction 
        #galcomps.AbsMag[maskmag] = CompCorMag - DistMod # No K correction applied
        #No K correction applied:
        galcomps.AbsMag[maskmag] = CompCorMag - dataned.DistMod2 # AbsMag using distance modulus independent of z 

        if ellconf.band in SunMag:
            MSun = SunMag[ellconf.band]

            #Lum = 10**((MSun - AbsMag)/2.5) Luminosity will  now be computed using AbsMag2
            datatidal.Lum = 10**((MSun - datatidal.AbsMag2)/2.5)
            # per component 
            galcomps.Lum[maskmag]= 10**((MSun - galcomps.AbsMag[maskmag])/2.5)

        else:
            print("Absolute Magnitude for Band {} was not found. Check filter name ".format(ellconf.band))
            print("Luminosity will not be computed.")

            datatidal.Lum = 0

        #print("Magnitud Absoluta",AbsMag)
        #print("Magnitud Absoluta using Distance Modulus independen of z ",AbsMag2)
        #print("check references in ",ellconf.namened)

    else:

        datatidal.AbsMag = 99
        datatidal.AbsMag2 = 99
        datatidal.Lum = 0


    if maskgalax.any():

        #if (ellconf.flagweb and ellconf.flagobj and not(ellconf.flagned)):

        galcomps.Rad50kpc[maskgalax] = galcomps.Rad50[maskgalax]*galhead.scale*dataned.Scalekpc

        galcomps.mme[maskgalax] = galcomps.mme[maskgalax] - dataned.GalExt - dataned.SbDim
        galcomps.me[maskgalax] = galcomps.me[maskgalax] - dataned.GalExt - dataned.SbDim
 
        #else:
        #    print("mean surface brightness at Re is not corrected for galactic extintion nor surface brightness dimming ")





    ################  INFORMATION CRITERIA ####################


    freepar = numParFree(galcomps) 
    npix = datatidal.ndof + freepar

    #    ;  AKAIKE INFORMATION CRITERION
    #    ; AIC = chi^2 + 2k

    datatidal.AICrit = datatidal.objchinu*datatidal.ndof + 2*freepar
    

    #    ; BAYESIAN INFORMATION CRITERION
    #    ; BIC = chi^2 + k * ln(n)

    datatidal.BICrit = datatidal.objchinu*datatidal.ndof + freepar*np.log(npix)


    #    ; BAYESIAN INFORMATION CRITERION limited by  
    # number of elements of resolution 

    #    ; BIC = chi^2 * APSF + k * ln(n/APSF)

    APSF = np.pi*ellconf.fwhm**2

    datatidal.BICres = datatidal.objchinu*datatidal.ndof + freepar*np.log(npix/APSF)

    ## for output only: 
    stidxg = np.argsort(sectgalax.radius)

    mgerad = sectgalax.radius[stidxg]

    #aell = mgerad.max() 
    #bell = mgerad.max()*ellconf.qarg



    #prints photometric variables into a file:
    printPhot(ellconf, galhead, galcomps, dataned, datatidal, sectgalax)

    passPhotVar(photapi, dataned, datatidal, galcomps)


    return photapi


def passPhotVar(photapi, dataned, datatidal, galcomps):

    photapi.aell = datatidal.aell
    photapi.bell = datatidal.bell 

    photapi.GalExt = dataned.GalExt 
    photapi.DistMod = dataned.DistMod
    photapi.DistMod2 = dataned.DistMod2
    photapi.Scalekpc = dataned.Scalekpc
    photapi.SbDim = dataned.SbDim 

    photapi.magalaper = datatidal.magalaper
    photapi.magmodaper=datatidal.magmodaper 
    photapi.totFlux = datatidal.totFlux 
    photapi.totMag = datatidal.totMag 
    photapi.BulgeToTotal = datatidal.BulgeToTotal
    photapi.tidal = datatidal.tidal 
    photapi.objchinu = datatidal.objchinu
    photapi.bump = datatidal.bump 
    photapi.snr = datatidal.snr
    photapi.stdsnr = datatidal.stdsnr 
    photapi.totsnr = datatidal.totsnr 
    photapi.rss   = datatidal.rss 
    photapi.ndof   = datatidal.ndof 
    photapi.AbsMag = datatidal.AbsMag 
    photapi.AbsMag2 = datatidal.AbsMag2 
    photapi.Lum = datatidal.Lum 
    photapi.AICrit   = datatidal.AICrit 
    photapi.BICrit = datatidal.BICrit
    photapi.BICres = datatidal.BICres

    photapi.EffRad = datatidal.EffRad
    photapi.EffRad9 = datatidal.EffRad9
    photapi.EffRad3 = datatidal.EffRad3

    photapi.KronRad = galcomps.KronRad.copy()
    photapi.PetRad = galcomps.PetRad.copy()





