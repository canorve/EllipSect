
from ellipsect.lib.libs import *

from ellipsect import *



def Tidal(params, galpar, galcomps, sectgalax, rmin):
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



    ab=galpar.q

    ell=1-ab

    aell = mgerad.max() 

    bell = mgerad.max() * ab

    #changing to arc sec
    aellarc=aell*galpar.scale

    #print("major axis, minor axis (pix) ",aell,bell)

    NCol=len(galpar.img[0])
    NRow=len(galpar.img)

    #print("max size ",NCol,NRow)

    #Obj.Angle = Obj.Theta - 90

    #angle computed from y-axis to x-axis  
    Theta=galpar.ang + 90

    if params.flagmodel == False:

        (xlo, xhi, ylo, yhi) = GetSize(galpar.xc, galpar.yc, aell, Theta, ab, NCol, NRow)

        #print("size box ",xmin, xmax, ymin, ymax)


        xser=galpar.xc
        yser=galpar.yc
    else:
        (xlo, xhi, ylo, yhi) = GetSize(galpar.inxc, galpar.inyc, aell, Theta, ab, NCol, NRow)

        xser=galpar.inxc
        yser=galpar.inyc




    imgal = galpar.img


    immodel = galpar.model


    imres = galpar.imres


    immask = galpar.mask


    hdu = fits.open(params.namesig)
    header=hdu[0].header  
    galpar.sigma=hdu[0].data
    #hdu.close()

    if params.flagmodel == False:
        imsigma = galpar.sigma.astype(float)
    else:
        imsigma = np.sqrt(np.abs(galpar.img)) #wrong but approx.
        masksigma = imsigma <= 0
        imsigma[masksigma] = 1 # avoiding division by zero
        print("I can't compute SNR. ")
        print("All SNR quantities will be wrong. ")


    # creates a new image for snr 
    #NCol=len(galpar.img[0])
    #NRow=len(galpar.img)
    #MakeImage(params.namesnr, NCol, NRow):


    header['TypeIMG'] = ('SNR', 'Signal to Noise Ratio image')
    hdu[0].header  =header
    galpar.imsnr=imgal/imsigma
    
    hdu[0].data = galpar.imsnr

    if params.flagsnr:
        hdu.writeto(params.namesnr, overwrite=True)
        print("SNR image created.. ",params.namesnr)
    hdu.close()


    #if galpar.tempmask!=None:
    imell=immask.copy()
    #else:
    #    imell=imgal.copy()


    imell.fill(False)


        #    for objchinu, Tidal and SNR
        #    maskm = dat[ylo - 1:yhi, xlo - 1:xhi] == num  # big image coordinates
        #maskm =immask[ylo - 1:yhi, xlo - 1:xhi] == False  # big image coordinates
    #if galpar.tempmask!=None:
    maskm =immask == False  # big image coordinates
    #else:
    #    maskm =imell == False  # big image coordinates


        #maskm =immask == False  # big image coordinates

        #   mask including rmin for Bumpiness only
        #    maskbum = dat[ylo - 1:yhi, xlo - 1:xhi] == num  # big image coordinates
        #maskbum = immask[ylo - 1:yhi, xlo - 1:xhi] == False  # big image coordinates
    #if galpar.tempmask!=None:
    maskbum = immask == False  # big image coordinates
    #else:
    #    maskbum = imell == False  # big image coordinates

    #############

    theta = 0

    ypos, xpos = np.mgrid[ylo - 1:yhi, xlo - 1:xhi]

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
    
    maskbum[ylo - 1:yhi, xlo - 1:xhi][mask] = False

    ## identifying area to compute photometry: 
    imell = ExtractEllip(imell, True, xser, yser, aell, Theta, ell, xlo, xhi, ylo, yhi)
    ## xlo, xhi, ylo, yhi makes sure that ellipse will not be outside of this range


    #maskm=maskm*imell
    #maskbum=maskbum*imell

    #maskm=maskm[ylo - 1:yhi, xlo - 1:xhi]*imell
    maskm=maskm*imell

    #maskbum=maskbum[ylo - 1:yhi, xlo - 1:xhi]*imell
    maskbum=maskbum*imell


    #if galpar.tempmask!=None:
    hdu = fits.open(galpar.tempmask)
    #else:
    #    hdu = fits.open(params.namesig)

    header=hdu[0].header  

    header['TypeIMG'] = ('Check', 'Use this image to check the area where photometry was computed')
    hdu[0].header  =header


    hdu[0].data = (~maskm).astype("int")*100
    hdu.writeto(params.namecheck, overwrite=True)
    hdu.close()




    if maskbum.any():

    # for Bumpiness

        #galfluxbum  = imgal[ylo - 1:yhi, xlo - 1:xhi][maskbum]
        #modfluxbum  = immodel[ylo - 1:yhi, xlo - 1:xhi][maskbum]
        #sigfluxbum  = imsigma[ylo - 1:yhi, xlo - 1:xhi][maskbum]
        galfluxbum  = imgal[maskbum]
        modfluxbum  = immodel[maskbum]
        sigfluxbum  = imsigma[maskbum]

    ####

    # for Tidal, SNR, and objchinu

    #        sumflux  = np.sum(imgal[ylo - 1:yhi, xlo - 1:xhi][maskm] - sky)
        #sumflux  = np.sum(imgal[ylo - 1:yhi, xlo - 1:xhi][maskm])
        #sumsig   = np.sum(imsigma[ylo - 1:yhi, xlo - 1:xhi][maskm])

        #galflux  = imgal[ylo - 1:yhi, xlo - 1:xhi][maskm]
        #modflux  = immodel[ylo - 1:yhi, xlo - 1:xhi][maskm]

        #sigflux  = imsigma[ylo - 1:yhi, xlo - 1:xhi][maskm]


        sumflux  = np.sum(imgal[maskm])
        sumfluxmod  = np.sum(immodel[maskm])
        sumsig   = np.sum(imsigma[maskm])


        galflux  = imgal[maskm]
        modflux  = immodel[maskm]

        sigflux  = imsigma[maskm]




        resflux = (galflux - modflux)**2

        #rss2=(resflux).sum()

        #print("Residual sum squares ",rss2)
        
        #  local chinu

        varchi = sigflux**2
        chinu  = np.sum(resflux/varchi)
        

        #pixcountchi = np.size(immodel[ylo - 1:yhi, xlo - 1:xhi][maskm])
        pixcountchi = np.size(immodel[maskm])


        if(pixcountchi > 11):

            ndof=pixcountchi - int(galcomps.freepar.sum())
            objchinu= chinu / ndof
        else:
            objchinu=-1
            ndof=-1


        # snr
        #if(np.size(galpar.imsnr[ylo - 1:yhi, xlo - 1:xhi][maskm]) > 0):

         #   totsnr=galpar.imsnr[ylo - 1:yhi, xlo - 1:xhi][maskm].sum()

          #  snr=galpar.imsnr[ylo - 1:yhi, xlo - 1:xhi][maskm].mean()
           # stdsnr=galpar.imsnr[ylo - 1:yhi, xlo - 1:xhi][maskm].std()

        if(np.size(galpar.imsnr[maskm]) > 0):

            totsnr=galpar.imsnr[maskm].sum()

            snr=galpar.imsnr[maskm].mean()
            stdsnr=galpar.imsnr[maskm].std()


        else:
            print("I can't compute SNR")
            snr=-1
            stdsnr=0
            totsnr=0
        # Tidal parameter

        tgal = np.abs((galflux)/(modflux) - 1)
        sumtidal=np.sum(tgal)
        #pixcountid=np.size(immodel[ylo - 1:yhi, xlo - 1:xhi][maskm])
        pixcountid=np.size(immodel[maskm])


        if pixcountid > 0:
            tidal = (sumtidal / pixcountid)
        else:
            tidal=-1

        # bumpiness

        resbump  = (galfluxbum - modfluxbum)**2
        # sflux    = np.sum(np.abs(modfluxbum - sky))
        sflux    = np.sum(np.abs(modfluxbum))
        varres   = sigfluxbum * sigfluxbum
        numbump  = np.sum(resbump - varres)

        #pixcountbum=np.size(immodel[ylo - 1:yhi, xlo - 1:xhi][maskbum])
        pixcountbum=np.size(immodel[maskbum])


        # Bumpiness
        if pixcountbum > 0:

            meansflux = sflux / pixcountbum
            meanres   = numbump / pixcountbum

            if (meanres < 0):
                meanres=0

            bump      = (np.sqrt(meanres)) / meansflux

        else:
            bump=-1


    # computing RSS: 
    #rss=(imres[ylo - 1:yhi, xlo - 1:xhi][maskm]**2).sum()
    rss=(imres[maskm]**2).sum()
    


    #computing magnitud for galaxy and model using aperture. 

    magalaper= galpar.mgzpt - 2.5*np.log10(sumflux/galpar.exptime)  #+ 0.1
    magmodaper = galpar.mgzpt - 2.5*np.log10(sumfluxmod/galpar.exptime) # + 0.1


  
    return (tidal,objchinu,bump,snr,stdsnr,totsnr,rss,ndof,magalaper,magmodaper)


#phot/tidal.py

def ExtractEllip(imagemat, idn, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
    "This subroutine extract an ellipse within an box to compute photometric parameters "
    "It returns the area of ellipse with idn values. "


    xmin = np.int(xmin)
    xmax = np.int(xmax)
    ymin = np.int(ymin)
    ymax = np.int(ymax)

    q = (1 - ell)
    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    ypos, xpos = np.mgrid[ymin - 1:ymax, xmin - 1:xmax]

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

    mask = dist < dell
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

    xmin = x - np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    xmax = x + np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    ymin = y - np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    ymax = y + np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    mask = xmin < 1
    if mask.any():
        if isinstance(xmin,np.ndarray):
            xmin[mask] = 1
        else:
            xmin = 1

    mask = xmax > ncol

    if mask.any():
        if isinstance(xmax,np.ndarray):
            xmax[mask] = ncol
        else:
            xmax = ncol

    mask = ymin < 1
    if mask.any():
        if isinstance(ymin,np.ndarray):
            ymin[mask] = 1
        else:
            ymin = 1

    mask = ymax > nrow
    if mask.any():
        if isinstance(ymax,np.ndarray):
            ymax[mask] = nrow
        else:
            ymax = nrow


    return (int(xmin), int(xmax), int(ymin), int(ymax))


