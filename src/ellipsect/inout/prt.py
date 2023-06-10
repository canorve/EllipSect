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

from ellipsect import *

import ellipsect

from ellipsect.inout.galfit  import numParFree

def PrintEllFilesGax(ellconf, galhead, xradq,ysbq,ysberrq,xradm,ysbm,ysberrm):
    "print surface brightness of galaxy and model to file"

    # output for galaxy
    filegal=ellconf.sboutput+".gal.txt"
    OUTFH = open (filegal,"w")

    lineout= "#        sectors_photometry used with q={} and pa={} (same as GALFIT) \n".format(ellconf.qarg, ellconf.parg)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {} \n".format(galhead.outimage, galhead.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galhead.exptime,galhead.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(ellconf.xc, ellconf.yc, ellconf.skylevel)
    OUTFH.write(lineout)

    lineout= "#            Galaxy                                   \n"
    OUTFH.write(lineout)

    lineout= "#     rad      SB        SBerr       \n"
    OUTFH.write(lineout)

    lineout= "# (arcsec) (mag/arcsec) (error)   \n"
    OUTFH.write(lineout)

    for idx, item in reversed(list(enumerate(xradq))):
        lineout= "{0:.3f} {1:.3f} {2:.3f} \n".format(xradq[idx],ysbq[idx],ysberrq[idx])
        OUTFH.write(lineout)

    OUTFH.close()

    # output for model 
    filemodel=ellconf.sboutput+".mod.txt"
    OUTFH = open (filemodel,"w")

    lineout= "#        sectors_photometry used with q={} and pa={} (same as GALFIT) \n".format(ellconf.qarg, ellconf.parg)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {} \n".format(galhead.outimage,galhead.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galhead.exptime,galhead.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(ellconf.xc, ellconf.yc, ellconf.skylevel)
    OUTFH.write(lineout)

    lineout= "#           Surface Brightness   Model                \n"
    OUTFH.write(lineout)

    lineout= "#    rad       SB        SBerr \n"
    OUTFH.write(lineout)

    lineout= "# (arcsec) (mag/arcsec) (error)  \n"
    OUTFH.write(lineout)

    for idx, item in reversed(list(enumerate(xradm))):
            
        lineout= "{0:.3f} {1:.3f} {2:.3f} \n".format(xradm[idx],ysbm[idx],ysberrm[idx])

        OUTFH.write(lineout)

    OUTFH.close()


    runcmd = "mv  {}  sbfiles/{}".format(filegal,filegal)
    errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)

    runcmd = "mv  {}  sbfiles/{}".format(filemodel,filemodel)
    errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)

#io/print.py
def PrintEllFilesComps(ellconf,galhead,namecomp,ncomp,xradq,ysbq,ysberrq):
    "Print surface brigthness of components to file "
    #subcomponent model 

    filesub = ellconf.sboutput+".comp-"+ncomp+".txt"
    OUTFH = open (filesub,"w")

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(ellconf.qarg,ellconf.parg)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {}  \n".format(galhead.outimage,galhead.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galhead.exptime,galhead.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(ellconf.xc, ellconf.yc, ellconf.skylevel)
    OUTFH.write(lineout)

    lineout= "#        Model  {}   component {}            \n".format(namecomp,ncomp)
    OUTFH.write(lineout)

    lineout= "#     rad      SB   SBerror  \n"
    OUTFH.write(lineout)

    lineout= "#   (arcsec) (mag/arcsec) (error) \n"
    OUTFH.write(lineout)

    for idx, item in reversed(list(enumerate(xradq))):

        lineout= "{0:.3f} {1:.3f} {2:.3f}  \n".format(xradq[idx],ysbq[idx],ysberrq[idx])

        OUTFH.write(lineout)

    OUTFH.close()

    runcmd = "mv  {}  sbfiles/{}".format(filesub,filesub)
    errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)


#io/print.py
def PrintFilesGax(ellconf,galhead,rtxtang,r,mgesb,angal,r2,mgemodsb,angmod):
    "Print surface parameters of galaxy and model to outfile "

    # galaxy
    filegalax=ellconf.sboutput+"-"+str(rtxtang)+".gal.txt"
    OUTFH = open (filegalax,"w")

    lineout= "# Surface brigthness profiles measured in ang = {} from major axis (0 deg)\n".format(rtxtang)
    OUTFH.write(lineout)

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(ellconf.qarg, ellconf.parg)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {} \n".format(galhead.outimage,galhead.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel]\n".format(galhead.exptime,galhead.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(ellconf.xc, ellconf.yc, ellconf.skylevel)
    OUTFH.write(lineout)

    lineout= "#            Galaxy                                   \n"
    OUTFH.write(lineout)

    lineout= "#     rad      SB              \n"
    OUTFH.write(lineout)

    lineout= "# (arcsec) (mag/arcsec)    \n"
    OUTFH.write(lineout)

    for idx, item in enumerate(r):

        lineout= "{0:.3f} {1:.3f} \n".format(r[idx],mgesb[angal][idx])

        OUTFH.write(lineout)

    OUTFH.close()

    #model
    filemodel=ellconf.sboutput+"-"+str(rtxtang)+".mod.txt"
    OUTFH = open (filemodel,"w")

    lineout= "# Surface brigthness profiles measured in ang = {} from major axis (0 deg)\n".format(rtxtang)
    OUTFH.write(lineout)

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(ellconf.qarg, ellconf.parg)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {} \n".format(galhead.outimage,galhead.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galhead.exptime,galhead.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(ellconf.xc, ellconf.yc, ellconf.skylevel)
    OUTFH.write(lineout)

    lineout= "#        Model                \n"
    OUTFH.write(lineout)

    lineout= "#     rad      SB     \n"
    OUTFH.write(lineout)

    lineout= "#   (arcsec) (mag/arcsec) \n"
    OUTFH.write(lineout)

    for idx, item in enumerate(r2):

        lineout= "{0:.3f} {1:.3f}  \n".format(r2[idx],mgemodsb[angmod][idx])

        OUTFH.write(lineout)

    OUTFH.close()

    runcmd = "mv  {}  sbfiles/{}".format(filegalax,filegalax)
    errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)

    runcmd = "mv  {}  sbfiles/{}".format(filemodel,filemodel)
    errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)


#io/print.py
def PrintFilesComps(ellconf,galhead,galcomps,rtxtang,ncomp,diffangle,rtemp,mgesbsub,ii,angtemp):

    #subcomponent model 

    filesub=ellconf.sboutput+"-"+str(rtxtang)+".comp-"+ncomp+".txt"
    OUTFH = open (filesub,"w")

    lineout= "# Surface brigthness profiles measured in ang = {} from major axis (0 deg)\n".format(rtxtang)
    OUTFH.write(lineout)

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(galcomps.AxRat[ii],90-galcomps.PosAng[ii])
    OUTFH.write(lineout)

    lineout= "# In the multiplot, there is a difference of {:.3f} for \n".format(diffangle)
    OUTFH.write(lineout)

    lineout= "# the one indicated in the top right corner.\n"
    OUTFH.write(lineout)

    lineout= "# The above is due to differences in the sectors_photometry.\n"
    OUTFH.write(lineout)

    lineout= "# for the galaxy and individual components.\n"
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {}  \n".format(galhead.outimage,galhead.mgzpt)
    OUTFH.write(lineout)

    lineout= "# exptime = {} plate scale = {} [arcsec per pixel] \n".format(galhead.exptime,galhead.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(ellconf.xc, ellconf.yc, ellconf.skylevel)
    OUTFH.write(lineout)

    lineout= "#        Model {} component {}            \n".format(galcomps.NameComp[ii],ncomp)
    OUTFH.write(lineout)

    lineout= "#     rad      SB     \n"
    OUTFH.write(lineout)

    lineout= "#   (arcsec) (mag/arcsec) \n"
    OUTFH.write(lineout)

    for idx, item in enumerate(rtemp):

        lineout= "{0:.3f} {1:.3f}  \n".format(rtemp[idx], mgesbsub[ii][angtemp][idx])

        OUTFH.write(lineout)

    OUTFH.close()

    runcmd = "mv  {}  sbfiles/{}".format(filesub,filesub)
    errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)


def printWelcome():


    print("Ellipsect: A surface brightness analysis tool for GALFIT output ")

    print("Version:",ellipsect.__version__)

    url = "https://github.com/canorve/EllipSect"
    
    print("webpage: "+url+"\n")


    print("check 'ellipsect -h' for full options \n")

    #print("EllipSect will create SB plots based on sectors ")
    #print("of an ellipse with the following properties: \n")


def printEllinfo(ellconf, galhead):



    print("center is at xc, yc = {}, {} ".format(ellconf.xc, ellconf.yc))

    str = "q = {} ".format(ellconf.qarg)
    print(str)

    str = "pa = {} ".format(ellconf.parg)
    print(str)


    ##
    str = "number of sectors = {}  ".format(ellconf.sectors)
    print(str)

    print("\nother parameters: \n")

    str = "Mag zeropoint = {} ".format(galhead.mgzpt)
    print(str)

    str = "Plate Scale = {} ".format(galhead.scale)
    print(str)

    str = "sky = {} ".format(ellconf.skylevel)
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


    print("In multi-plot, every color represents the same as the ones in the single-plot's legend")




def printPhot(ellconf, galhead, galcomps, dataned, datatidal, sectgalax):


    #redefining the mask again
    masknuker = (galcomps.NameComp == "nuker") & (galcomps.Active == True)

    maskmag = ((galcomps.NameComp != "ferrer") & (galcomps.NameComp != "nuker") 
                & (galcomps.NameComp != "edgedisk") & (galcomps.NameComp != "king") 
                & (galcomps.Active == True))





    #note separate this in three functions one for the header
    # other for the output variables for file
    # and other for the output component variables for the class

    #separate this in another function
    #######  file output:  ######

    print("Creating output photometry file: ",ellconf.output)

    OUTPHOT = open (ellconf.output,"w")

    lineout= "#   Output photometry for {} file \n".format(galhead.outimage)
    OUTPHOT.write(lineout)

    lineout= "#\n"
    OUTPHOT.write(lineout)

    lineout = "# Photometry for object: {} (use -object option to change it) \n".format(ellconf.objname)
    OUTPHOT.write(lineout)

    lineout = "# All photometric quantities are computed for filter {} (use -filter option to change it) \n".format(ellconf.band)
    OUTPHOT.write(lineout)



    lineout = '# Cosmology: hconst = {}, omega m = {}, omega lambda= {}'.format(ellconf.hconst,ellconf.omegam,ellconf.omegav) 
    OUTPHOT.write(lineout)


    lineout = "# Magnitudes are not corrected by K-Correction \n"
    OUTPHOT.write(lineout)

    lineout= "# Total magnitude and other variables are NOT computed for: \n"
    OUTPHOT.write(lineout)

    lineout= "# ferrer, nuker, edgedisk and king components.  \n"
    OUTPHOT.write(lineout)

    if masknuker.any():
 
        lineout= "# For Nuker component, Rad90 is the gamma radius   \n"
        OUTPHOT.write(lineout)



    lineout = "# Some of the photometric variables mentioned below were computed within an ellipse containing 90% of total light \n"
    OUTPHOT.write(lineout)

    lineout = "# This ellipse has axis a = {:.2f} and b = {:.2f} centered at xc, yc \n".format(datatidal.aell, datatidal.bell)
    OUTPHOT.write(lineout)

    lineout = "# To see this ellipse check the file {} \n".format(ellconf.namecheck)
    OUTPHOT.write(lineout)


    lineout= "#\n"
    OUTPHOT.write(lineout)

    lineout= "# sectors_photometry was used with q = {:.2f}, pa = {:.2f} and minlevel = {} \n".format(ellconf.qarg, 
                                                                                            ellconf.parg, ellconf.minlevel)
    OUTPHOT.write(lineout)

    lineout= "# OutImage = {}  MgZpt = {}  \n".format(galhead.outimage,galhead.mgzpt)
    OUTPHOT.write(lineout)

    lineout= "# exptime = {}  plate scale = {} ''/pix \n".format(galhead.exptime,galhead.scale)
    OUTPHOT.write(lineout)

    lineout= "# xc = {:.2f}  yc = {:.2f}  sky = {:.2f} \n".format(ellconf.xc, ellconf.yc, ellconf.skylevel)
    OUTPHOT.write(lineout)
    
    lineout = "# Gal. Extinction = {}; Distance Mod. = {}; Distance Mod. (z independent) = {} \n".format(
                    dataned.GalExt, dataned.DistMod, dataned.DistMod2)
    OUTPHOT.write(lineout)

    lineout = "# cosmology corrected scale = {} kpc/arcsec; Surface brightness dimming (mag) = {} \n".format(
                        dataned.Scalekpc, dataned.SbDim)
    OUTPHOT.write(lineout)

    lineout = "# Check references in the NED file: {} \n\n".format(ellconf.namened)
    OUTPHOT.write(lineout)

    #note this below goes in another function

    lineout = "total apparent mag of the galaxy (using ellipse aperture) = {:.3f}  \n".format(datatidal.magalaper)
    OUTPHOT.write(lineout)

    lineout = "total apparent mag of the model (using ellipse aperture) = {:.3f}  \n".format(datatidal.magmodaper)
    OUTPHOT.write(lineout)


    if maskmag.any():
        #Flux=10**((galhead.mgzpt -  galcomps.Mag[maskmag])/2.5)
        #totFlux=Flux.sum()
        #datatidal.totMag=-2.5*np.log10(totFlux) + galhead.mgzpt

        lineout = "total apparent mag from model parameters (without corrections) = {:.3f}  \n".format(datatidal.totMag)
        OUTPHOT.write(lineout)

        lineout = "total flux (without corrections) = {:.3f} \n\n".format(datatidal.totFlux)
        OUTPHOT.write(lineout)
    else:
        lineout="Total magnitude can not be computed with the actual galfit functions \n\n"
        OUTPHOT.write(lineout)


    lineout="Bulge To Total Ratio = {:.3f} \n".format(datatidal.BulgeToTotal)
    OUTPHOT.write(lineout)

    lineout="Tidal = {:.3f} \n".format(datatidal.tidal)
    OUTPHOT.write(lineout)

    lineout="Local Chinu = {:.3f} \n".format(datatidal.objchinu)
    OUTPHOT.write(lineout)

    lineout="Bumpiness = {:.3f} \n".format(datatidal.bump)
    OUTPHOT.write(lineout)

    lineout = "mean SNR = {:.3f} \n".format(datatidal.snr)
    OUTPHOT.write(lineout)

    lineout = "std SNR = {:.3f} \n".format(datatidal.stdsnr)
    OUTPHOT.write(lineout)

    lineout = "total SNR sum = {:.3f} \n".format(datatidal.totsnr)  
    OUTPHOT.write(lineout)

    lineout = "Residual sum of squares (RSS) = {:.3f}  \n".format(datatidal.rss)  
    OUTPHOT.write(lineout)

    lineout= "degrees of freedom = {} \n".format(datatidal.ndof)  
    OUTPHOT.write(lineout)

    maskfreegal = galcomps.Active == True

    totfreepar = numParFree(galcomps) 

    lineout="Number of free params = {} \n\n".format(totfreepar)
    OUTPHOT.write(lineout)
   

    if maskmag.any():


        if (ellconf.flagweb and ellconf.flagobj and not(ellconf.flagned)):

            #CorMag = totMag - GalExt # corrected by galactic extinction 
            
            #AbsMag=CorMag - DistMod # No K correction applied

            #AbsMag2=CorMag - DistMod2 # No K correction applied

            lineout="Absolute Mag = {:.3f} \n".format(datatidal.AbsMag)
            OUTPHOT.write(lineout)

            lineout="Absolute Mag using Dist Mod independent of z = {:.3f} \n".format(datatidal.AbsMag2)
            OUTPHOT.write(lineout)

            lineout="Luminosity = {:.3f} (10^10 solar lum) using Dist Mod independent of z  \n\n".format(datatidal.Lum/1e10)
            OUTPHOT.write(lineout)



    lineout= "Akaike Information Criterion = {:.3f} \n".format(datatidal.AICrit)  
    OUTPHOT.write(lineout)

    lineout = "Bayesian Information Criterion = {:.3f} \n".format(datatidal.BICrit)  
    OUTPHOT.write(lineout)

    lineout = "Bayesian Information Criterion using nres = n / Area_psf = {:.3f} \n".format(datatidal.BICres)  
    OUTPHOT.write(lineout)

    lineout= "\n"
    OUTPHOT.write(lineout)



    lineout= "The following Radius were computed using all the " \
            + "Sersic (and related) components:\n"
    OUTPHOT.write(lineout)


    lineout= "Radius at 30% of the total galaxy light: {:.3f} \n".format(datatidal.EffRad3)  
    OUTPHOT.write(lineout)

    lineout= "Effective Radius of the total galaxy light: {:.3f} \n".format(datatidal.EffRad)  
    OUTPHOT.write(lineout)

    lineout= "Radius at 90% of the total galaxy light: {:.3f} \n".format(datatidal.EffRad9)  
    OUTPHOT.write(lineout)




    lineout= "\n"
    OUTPHOT.write(lineout)

    ############################

    if ellconf.flagradsky:

        mean = ellconf.gradskymean
        std = ellconf.gradskystd
        median = ellconf.gradskymed


        lineout="computed sky with the gradient method:  mean = {:.2f} , std = {:.2f}, median = {} \n".format(mean,std,median)

        OUTPHOT.write(lineout)
 
    if ellconf.flagrandboxsky:

        mean = ellconf.randskymean
        std = ellconf.randskystd
        median = ellconf.randskymed

        lineout="computed sky with the random box method:  mean = {:.2f} , std = {:.2f}, median = {} \n".format(mean,std,median)
        OUTPHOT.write(lineout)


    lineout = "\n"  
    OUTPHOT.write(lineout)

    #############################
    #note this goes in another functions

    lineout = "#########################################\n"  
    OUTPHOT.write(lineout)

    lineout = "# Photometric properties per component: #\n"  
    OUTPHOT.write(lineout)

    lineout = "#########################################\n"  
    OUTPHOT.write(lineout)

    lineout= "\n"
    OUTPHOT.write(lineout)


    lineout = "#   N      Component   FractionLight    me        <me>        AbsMag      Luminosity      Rad90       Re      KronRadius   PetroRad \n"  
    OUTPHOT.write(lineout)

    lineout = "#                                    (mag/'')   (mag/'')      (mag)    (10^10 SolarLum)   (pix)      (kpc)      (pix)       (Re)    \n"  
    OUTPHOT.write(lineout)


    for idx, item in enumerate(galcomps.N) :
        lineout= "    {0:^2} {1:^17} {2:^10.3f} {3:^10.3f} {4:^10.3f} {5:^14.3f} {6:^14.3f} {7:^10.3f} {8:^10.3f} {9:^10.3f}  {10:^10.3f}  \n".format(galcomps.N[idx],galcomps.NameComp[idx],galcomps.PerLight[idx],galcomps.me[idx],galcomps.mme[idx],galcomps.AbsMag[idx],galcomps.Lum[idx]/1e10,galcomps.Rad90[idx],galcomps.Rad50kpc[idx],galcomps.KronRad[idx],galcomps.PetRad[idx])
        OUTPHOT.write(lineout)

    OUTPHOT.close()

 











