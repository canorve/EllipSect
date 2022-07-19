
from ellipsect.lib.libs import *

from ellipsect import *

import ellipsect


def PrintEllFilesGax(params,galpar,xradq,ysbq,ysberrq,xradm,ysbm,ysberrm):
    "print surface brightness of galaxy and model to file"

    # output for galaxy
    filegal=params.sboutput+".gal.txt"
    OUTFH = open (filegal,"w")

    lineout= "#        sectors_photometry used with q={} and pa={} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {} \n".format(galpar.outimage,galpar.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galpar.exptime,galpar.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
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
    filemodel=params.sboutput+".mod.txt"
    OUTFH = open (filemodel,"w")

    lineout= "#        sectors_photometry used with q={} and pa={} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {} \n".format(galpar.outimage,galpar.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galpar.exptime,galpar.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
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
def PrintEllFilesComps(params,galpar,namecomp,ncomp,xradq,ysbq,ysberrq):
    "Print surface brigthness of components to file "
    #subcomponent model 

    filesub = params.sboutput+".comp-"+ncomp+".txt"
    OUTFH = open (filesub,"w")

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {}  \n".format(galpar.outimage,galpar.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galpar.exptime,galpar.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
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
def PrintFilesGax(params,galpar,rtxtang,r,mgesb,angal,r2,mgemodsb,angmod):
    "Print surface parameters of galaxy and model to outfile "

    # galaxy
    filegalax=params.sboutput+"-"+str(rtxtang)+".gal.txt"
    OUTFH = open (filegalax,"w")

    lineout= "# Surface brigthness profiles measured in ang = {} from major axis (0 deg)\n".format(rtxtang)
    OUTFH.write(lineout)

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {} \n".format(galpar.outimage,galpar.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel]\n".format(galpar.exptime,galpar.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
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
    filemodel=params.sboutput+"-"+str(rtxtang)+".mod.txt"
    OUTFH = open (filemodel,"w")

    lineout= "# Surface brigthness profiles measured in ang = {} from major axis (0 deg)\n".format(rtxtang)
    OUTFH.write(lineout)

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {} \n".format(galpar.outimage,galpar.mgzpt)
    OUTFH.write(lineout)

    lineout= "#  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galpar.exptime,galpar.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
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
def PrintFilesComps(params,galpar,galcomps,rtxtang,ncomp,diffangle,rtemp,mgesbsub,ii,angtemp):

    #subcomponent model 

    filesub=params.sboutput+"-"+str(rtxtang)+".comp-"+ncomp+".txt"
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

    lineout= "#  OutImage = {}  magzpt = {}  \n".format(galpar.outimage,galpar.mgzpt)
    OUTFH.write(lineout)

    lineout= "# exptime = {} plate scale = {} [arcsec per pixel] \n".format(galpar.exptime,galpar.scale)
    OUTFH.write(lineout)


    lineout= "#  xc = {:.2f}  yc = {:.2f}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
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

    print("EllipSect will create SB plots based on sectors ")
    print("of an ellipse with the following properties: \n")







