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

from ellipsect.inout.galfit  import numParFree

def Tidal(datatidal, ellconf, dataimg, galhead, galcomps, sectgalax, rmin):
    "Computes Tidal  values as defined in Tal et al. 2009 AJ. It algo computes Bumpiness"
    "(Blakeslee 2006 ApJ) value defined between rmin and radius"

    # rmin = minimum radius to compute Tidal and Bumpiness (to avoid PSF Mismatch)

    #init values

    pixcount=0
    pixcountsnr=0
    pixcountchi=0
    sumtidal=0
    sumsig=0
    sigma=1
    flux=0
    sumflux=0
    meanflux=0
    resbump=0
    sumres=0
    sumchinu=0
    chinu=0
    varchi=1
    sflux=0
    meanres=0
    varres=0
    numbump=0
    ndof=0

    rss=0

    tidal=0
    objchinu=0
    bump=0
    snr=0

    totsnr=0
    stdsnr=0

    magalaper= 99
    magmodaper = 99


    ################

  
    stidxg = np.argsort(sectgalax.radius)

    mgerad=sectgalax.radius[stidxg]
    mgecount=sectgalax.counts[stidxg]
    mgeangle=sectgalax.angle[stidxg]
    mgeanrad=np.deg2rad(mgeangle)



    ab = ellconf.qarg
    ell = 1 - ab


    aell = datatidal.aell
    bell = datatidal.bell

    


    #changing to arc sec
    aellarc = aell*galhead.scale

    #print("major axis, minor axis (pix) ",aell,bell)

    NCol = len(dataimg.img[0])
    NRow = len(dataimg.img)

    #print("max size ",NCol,NRow)

    #Obj.Angle = Obj.Theta - 90

    #angle computed from y-axis to x-axis  
    Theta = ellconf.parg + 90

    if ellconf.flagmodel == False:

        (xlo, xhi, ylo, yhi) = GetSize(ellconf.xc, ellconf.yc, aell, Theta, ab, NCol, NRow)

        #print("size box ",xmin, xmax, ymin, ymax)


        xser = ellconf.xc
        yser = ellconf.yc
    else:
        #note: Ncol and Nrow are wrong in this situation because correspond to the large img
        (xlo, xhi, ylo, yhi) = GetSize(ellconf.inxc, ellconf.inyc, aell, Theta, ab, NCol, NRow)

        xser = ellconf.inxc
        yser = ellconf.inyc




    imgal = dataimg.img


    immodel = dataimg.model


    imres = dataimg.imres


    immask = dataimg.mask


    hdu = fits.open(ellconf.namesig)
    header = hdu[0].header  
    dataimg.sigma = hdu[0].data
    #hdu.close()

    if ellconf.flagmodel == False:
        imsigma = dataimg.sigma.astype(float)
    else:
        imsigma = np.sqrt(np.abs(dataimg.img)) #wrong but approx.
        masksigma = imsigma <= 0
        imsigma[masksigma] = 1 # avoiding division by zero
        print("I can't compute SNR. ")
        print("All SNR quantities will be wrong. ")


    # creates a new image for snr 
    #NCol=len(dataimg.img[0])
    #NRow=len(dataimg.img)
    #MakeImage(ellconf.namesnr, NCol, NRow):

    ##################################
    # creation of the  SNR image 
    ##################################

    header['TypeIMG'] = ('SNR', 'Signal to Noise Ratio image')
    hdu[0].header  =header


    if imgal.shape == imsigma.shape:

        dataimg.imsnr=imgal/imsigma

    else:
        print("WARNING: SNR can not be computed because"
                + "galaxy and sigma images have different shapes")
        dataimg.imsnr=np.ones(imgal.shape)
    
    hdu[0].data = dataimg.imsnr

    if ellconf.flagsnr:
        hdu.writeto(ellconf.namesnr, overwrite=True)
        print("SNR image created.. ",ellconf.namesnr)


    ##################################
    # creation of the Chi-square image 
    ##################################

    #if ellconf.flagchi:

    header['TypeIMG'] = ('CHI-SQR', 'Chi-Square image')
    hdu[0].header  = header


    maskchi = immask == True # big image coordinates


    if imgal.shape == imsigma.shape:
        dataimg.imchi = (imgal - immodel)**2 / imsigma**2
        if maskchi.any():
            dataimg.imchi[maskchi] = 0
    else:
        print("WARNING: Chi-square image can not be computed" +
                "because galaxy and sigma images have different shapes")
        dataimg.imchi = np.ones(imgal.shape)
    
    hdu[0].data = dataimg.imchi

    hdu.writeto(ellconf.namechi, overwrite=True)
    print("Chi square image created.. ", ellconf.namechi)


    ##################################
    # end of the Chi-square image 
    ##################################




    hdu.close()


    #if galhead.tempmask!=None:
    imell = immask.copy()
    #else:
    #    imell=imgal.copy()


    imell.fill(False)


    #if galhead.tempmask!=None:
    maskm = immask == False  # big image coordinates
    #else:
    #    maskm =imell == False  # big image coordinates


        #maskm =immask == False  # big image coordinates

        #   mask including rmin for Bumpiness only
    #if galhead.tempmask!=None:
    maskbum = immask == False  # big image coordinates
    #else:
    #    maskbum = imell == False  # big image coordinates

    #############

    theta = 0

    ypos, xpos = np.mgrid[ylo - 1: yhi + 1 , xlo - 1: xhi + 1 ]

    dx = xpos - xser
    dy = ypos - yser

    dist = np.sqrt(dx**2 + dy**2)

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / rmin, np.cos(landa) / rmin)

    xell = xser + rmin * np.cos(angle)
    yell = yser + rmin * np.sin(angle)

    dell = np.sqrt((xell - xser)**2 + (yell - yser)**2)

    mask = dist < dell

    #  correcting for rmin
    
    maskbum[ylo - 1: yhi + 1, xlo - 1: xhi +1][mask] = False


    ## identifying area to compute photometry: 
    imell = ExtractEllip(imell, True, xser, yser, aell, Theta, ell, xlo, xhi, ylo, yhi)
    ## xlo, xhi, ylo, yhi makes sure that ellipse will not be outside of this range


    # future coder: in case you doubt that the 
    #following two lines are incorrect: These are not. I double check
    maskm = maskm*imell 
    maskbum = maskbum*imell



    hdu = fits.open(galhead.tempmask)

    header = hdu[0].header  

    header['TypeIMG'] = ('Check', 'Use this image to check the' +  
                            'area where photometry was computed')
    hdu[0].header = header


    hdu[0].data = (~maskm).astype("int")*100
    hdu.writeto(ellconf.namecheck, overwrite=True)
    hdu.close()




    if maskbum.any():

        galfluxbum  = imgal[maskbum]
        modfluxbum  = immodel[maskbum]
        if imgal.shape == imsigma.shape:

            sigfluxbum  = imsigma[maskbum]
            sumsig   = np.sum(imsigma[maskm])
            sigflux  = imsigma[maskm]
        else:
            print("WARNING: Bumpiness and local chinu will have a wrong values because galaxy and sigma images have different shapes")
 
            sigfluxbum  = np.ones(imgal[maskbum].shape)
            sumsig = np.ones(imgal[maskm].shape)
            sigflux= np.ones(imgal[maskm].shape)


    ####

    # for Tidal, SNR, and objchinu

        sumflux  = np.sum(imgal[maskm])
        sumfluxmod  = np.sum(immodel[maskm])


        galflux  = imgal[maskm]
        modflux  = immodel[maskm]





        resflux = (galflux - modflux)**2

       
        #  local chinu

        varchi = sigflux**2
        chinu  = np.sum(resflux/varchi)
        

        pixcountchi = np.size(immodel[maskm])

        totfreepar = numParFree(galcomps) 

        if(pixcountchi > 11):

            ndof = pixcountchi - totfreepar 
            objchinu = chinu/ndof
        else:
            objchinu =- 1
            ndof =- 1


        if(np.size(dataimg.imsnr[maskm]) > 0):

            totsnr = dataimg.imsnr[maskm].sum()

            snr = dataimg.imsnr[maskm].mean()
            stdsnr = dataimg.imsnr[maskm].std()


        else:
            print("I can't compute SNR")
            snr=-1
            stdsnr=0
            totsnr=0
        # Tidal parameter

        tgal = np.abs((galflux)/(modflux) - 1)
        sumtidal = np.sum(tgal)

        pixcountid = np.size(immodel[maskm])


        if pixcountid > 0:
            tidal = (sumtidal / pixcountid)
        else:
            tidal=-1

        # bumpiness

        resbump  = (galfluxbum - modfluxbum)**2
        sflux    = np.sum(np.abs(modfluxbum))
        varres   = sigfluxbum * sigfluxbum
        numbump  = np.sum(resbump - varres)

        pixcountbum = np.size(immodel[maskbum])


        # Bumpiness
        if pixcountbum > 0:

            meansflux = sflux / pixcountbum
            meanres   = numbump / pixcountbum

            if (meanres < 0):
                meanres = 0

            bump      = (np.sqrt(meanres)) / meansflux

        else:
            bump=-1


        #ks statistics
        #ks,p = stats.kstest(galflux.flatten(),modflux.flatten())
        #print("kolmogorov ks, p: ",ks,p)

    # computing RSS: 
    rss = (imres[maskm]**2).sum()
    


    #computing magnitud for galaxy and model using aperture. 

    magalaper = galhead.mgzpt - 2.5*np.log10(sumflux/galhead.exptime)  #+ 0.1
    magmodaper = galhead.mgzpt - 2.5*np.log10(sumfluxmod/galhead.exptime) # + 0.1


    datatidal.tidal = tidal
    datatidal.objchinu = objchinu
    datatidal.bump = bump
    datatidal.snr = snr 
    datatidal.stdsnr = stdsnr
    datatidal.totsnr = totsnr
    datatidal.rss = rss
    datatidal.ndof = ndof
    datatidal.magalaper = magalaper
    datatidal.magmodaper = magmodaper 
  
    return datatidal 


#phot/tidal.py

def ExtractEllip(imagemat, idn, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
    "This subroutine extract an ellipse within an box to compute photometric parameters "
    "It returns the area of ellipse with idn values. "


    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    q = (1 - ell)
    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    ypos, xpos = np.mgrid[ymin - 1: ymax + 1, xmin - 1: xmax + 1]

    dx = xpos - x
    dy = ypos - y

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

    xell = x + R * np.cos(angle) * np.cos(theta) - bim * \
        np.sin(angle) * np.sin(theta)
    yell = y + R * np.cos(angle) * np.sin(theta) + bim * \
        np.sin(angle) * np.cos(theta)

    dell = np.sqrt((xell - x)**2 + (yell - y)**2)
    dist = np.sqrt(dx**2 + dy**2)

    mask = dist <= dell
    imagemat[ypos[mask], xpos[mask]] = idn

    return imagemat


#phot/tidal.py
def GetSize(x, y, R, theta, q, ncol, nrow):
    "this subroutine get the maximun"
    "and minimim pixels for Kron and sky ellipse"
    
    #theta is measured from x-axis    
    #q = (1 - ell)
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!


    # getting size

    constx =  np.sqrt((R**2)*(np.cos(theta))**2 + (bim**2)*(np.sin(theta))**2)
    consty =  np.sqrt((R**2)*(np.sin(theta))**2 + (bim**2)*(np.cos(theta))**2)



    xmin = x - constx                       
    xmax = x + constx
    ymin = y - consty
    ymax = y + consty
                    


    mask = xmin < 1
    if mask.any():
        if isinstance(xmin,np.ndarray):
            xmin[mask] = 1
        else:
            xmin = 1

    mask = xmax > ncol

    if mask.any():
        if isinstance(xmax,np.ndarray):
            xmax[mask] = ncol - 1 
        else:
            xmax = ncol - 1

    mask = ymin < 1
    if mask.any():
        if isinstance(ymin,np.ndarray):
            ymin[mask] = 1
        else:
            ymin = 1

    mask = ymax > nrow
    if mask.any():
        if isinstance(ymax,np.ndarray):
            ymax[mask] = nrow - 1
        else:
            ymax = nrow - 1


    return (round(xmin), round(xmax), round(ymin), round(ymax))


