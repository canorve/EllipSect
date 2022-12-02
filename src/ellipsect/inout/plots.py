
from ellipsect.lib.libs import *

from ellipsect import *



def PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,ellconf,scale):
    """  Produces final best-fitting plot  """

    # subplot for arc sec axis
    plt.close('all')


    #ULISES begin
    fig, (axsec,axred) = plt.subplots(2, sharex=True, sharey=False)
    gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
    gs.update(hspace=0.07)
    #ULISES end 


    if ellconf.flagranx == True:
        (xmin,xmax)=ellconf.ranx[0], ellconf.ranx[1]

    if ellconf.flagrany == True:
        (ymin,ymax)=ellconf.rany[0], ellconf.rany[1]


    minrad = np.min(xradq)
    maxrad = np.max(xradq)

    mincnt = np.min(ysbq)
    maxcnt = np.max(ysbq)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
    yran = mincnt * (maxcnt/mincnt)**np.array([+1.05, -0.05]) #inverted axis

   
    # ULISES begin
    axsec = plt.subplot(gs[0])

    axsec.set_ylabel("Surface Brightness (mag/'')")
    # ULISES end

    axsec.errorbar(xradq, ysbq,yerr=ysberrq,fmt='o-',capsize=2,color='red',markersize=0.7,label="galaxy",linewidth=2)


    if ellconf.flagalax == False:
        axsec.errorbar(xradm, ysbm,yerr=ysberrm,fmt='o-',capsize=2,color='blue',markersize=0.7,label="Model",linewidth=2)


    if ellconf.flagrany == True:
        axsec.set_ylim(ymax,ymin) #inverted
    else:
        axsec.set_ylim(yran)

    if ellconf.flaglogx == True:

        axsec.set_xscale("log")

        locmaj = LogLocator(base=10,numticks=12)
        axsec.xaxis.set_major_locator(locmaj)

        locmin = LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
        axsec.xaxis.set_minor_locator(locmin)
        axsec.xaxis.set_minor_formatter(NullFormatter())

    else:
        axsec.xaxis.set_minor_locator(AutoMinorLocator())
        axsec.xaxis.set_major_locator(AutoLocator())


    axsec.tick_params(which='both', width=2)
    axsec.tick_params(which='major', length=7)
    axsec.tick_params(which='minor', length=4, color='r')


    #begin psf fwhm 
    if ellconf.flagfwhm: 
        xpos = ellconf.fwhm*scale
        axsec.axvline(x=xpos,  linestyle='--', color='k', linewidth=2)
    # end 


    # ULISES begin
    axsec.axes.xaxis.set_ticklabels([])
    # ULISES end 


    axsec.yaxis.set_minor_locator(AutoMinorLocator())
    axsec.yaxis.set_major_locator(AutoLocator())


    # ULISES begin
    # Residual plot
    if len(ysbq) < len(ysbm):
        ysbm = ysbm[len(ysbm)-len(ysbq):]
    
    elif len(ysbq) > len(ysbm):
        ysbq = ysbq[len(ysbq)-len(ysbm):]
        ysberrq = ysberrq[len(ysberrq)-len(ysbq):]

    if len(ysbq) < len(ysberrm):
        ysberrm = ysberrm[len(ysberrm)-len(ysbq):]
    elif len(ysbq) > len(ysberrm):
        ysbq = ysbq[len(ysbq)-len(ysberrm):]        

    residual = ((ysbq-ysbm)/ysbq)*100 # (data-model)/data in percentage

    err = ((ysbm/ysbq**2)**2) * ysberrq**2 + ((1/ysbq)**2) * ysberrm**2 
    err = np.sqrt(err)*100
    axred = plt.subplot(gs[1])

    if len(xradq) != len(residual):
        axred.errorbar(xradm, residual, yerr=err,fmt='.',capsize=2,color='k')
    else:
        axred.errorbar(xradq, residual, yerr=err,fmt='.',capsize=2,color='k')

    axred.axhline(y=0,ls='dashed', color='k')
    axred.set_xlabel('Radius (arcsec)')
    axred.set_ylabel('Residual (%)')
    axred.set_ylim(-2,2)
    # ULISES end


    if ellconf.flaglogx == True:

        axred.set_xscale("log")

        locmaj = LogLocator(base=10,numticks=12)
        axred.xaxis.set_major_locator(locmaj)

        locmin = LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
        axred.xaxis.set_minor_locator(locmin)
        axred.xaxis.set_minor_formatter(NullFormatter())
    else:
        axred.xaxis.set_minor_locator(AutoMinorLocator())
        axred.xaxis.set_major_locator(AutoLocator())



    axred.tick_params(which='both', width=2)
    axred.tick_params(which='major', length=7)
    axred.tick_params(which='minor', length=4, color='r')

    if ellconf.flagranx == True:
        axsec.set_xlim(xmin,xmax)
        axred.set_xlim(xmin,xmax) #ulises plot
    else:
        axsec.set_xlim(xran)
        axred.set_xlim(xran) #ulises plot
 



    if ellconf.flagpix == True:

        axpix = axsec.twiny()

        axpix.set_xlabel("Pixels")
        x1, x2 = axsec.get_xlim()

        axpix.set_xlim(x1/scale, x2/scale)

        axpix.figure.canvas.draw()


        if ellconf.flaglogx == True:
            axpix.set_xscale("log")
            locmaj = LogLocator(base=10,numticks=12)
            axpix.xaxis.set_major_locator(locmaj)

            locmin = LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
            axpix.xaxis.set_minor_locator(locmin)
            axpix.xaxis.set_minor_formatter(NullFormatter())
        else:
            axpix.xaxis.set_minor_locator(AutoMinorLocator())
            axpix.xaxis.set_major_locator(AutoLocator())

        axpix.tick_params(which='both', width=2)
        axpix.tick_params(which='major', length=7)
        axpix.tick_params(which='minor', length=4, color='r')

        axret=axsec
    else:
        axret=axsec

    if ellconf.flagrid == True:
        # Customize the major grid
        axsec.grid(which='major', linestyle='-', linewidth='0.7', color='black')
        # Customize the minor grid
        axsec.grid(which='minor', linestyle=':', linewidth='0.5', color='black')


    #change the linewidth of the axis
    for axis in ['top','bottom','left','right']:
        axsec.spines[axis].set_linewidth(1.5)
        axred.spines[axis].set_linewidth(1.5)


    return xran,yran,axret

#io/plot.py
def PlotSub(xradq,ysbq,nsub,axsec,namec,colorval):
    """
    Produces subcomponent plot

    """

    substr=namec + " " + np.str_(nsub+1)

    axsec.plot(xradq, ysbq,'--',color=colorval,linewidth=1.7,markersize=0.7,label=substr)
    #axsec.plot(xradq, ysbq,'--',color=colorval,linewidth=1.5,markersize=0.7,label=substr)


def plotCube(ellconf, galhead, galcomps):

    ######################
    #shows the image cube#
    ######################{

    linewidth = 1.2

    if (ellconf.flagcomp):
        ell = Comp2Ellip(galhead, galcomps, ellconf.tot, linewidth)
    else:
        ell=[]


    ShowCube(galhead.outimage, namepng = ellconf.namecube, dpival 
             = ellconf.dpival, bri = ellconf.brightness, con = ellconf.contrast, 
             frac = ellconf.frac, fracmax = ellconf.fracmax,  cmap = ellconf.cmap, 
             ellipse = ell)


    if ellconf.dplot:
        plt.pause(1.5)
 
    plt.close()

    #}
    #####################


def Comp2Ellip(galhead, galcomps, N, lw=1):
    ''' converts galfit component parameter into an Ellipse object''' 


    ellipses = [] 

    #color value
    values = range(N)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    
    for idx, item in enumerate(galcomps.N):

        if galcomps.Activate[idx] == True:
            # correcting coordinates
            xc = galcomps.PosX[idx] - galhead.xmin + 1
            yc = galcomps.PosY[idx] - galhead.ymin + 1


            pa = galcomps.PosAng[idx] + 90

            w = galcomps.Rad[idx]
            h = galcomps.Rad[idx]*galcomps.AxRat[idx]


            colorVal = scalarMap.to_rgba(values[idx])


            ell=Ellipse((xc, yc), width = w, height = h, angle = pa,
                         edgecolor = colorVal,
                         facecolor = 'none',
                         linewidth = lw)

            ellipses.append(ell)


    return ellipses






class ShowCube:

    def __init__(self, cubeimg: str, namepng="cubeout.png", dpival=100, 
                bri = 33, con = 0.98, frac = 1, fracmax = 1, cmap='viridis', ellipse=[]):
        """
        This routine shows the GALFIT output cube image: galaxy, model and residual    
        """

        
        hdu = fits.open(cubeimg)
        data = (hdu[1].data.copy()).astype(float)
        model = (hdu[2].data.copy()).astype(float)
        residual = (hdu[3].data.copy()).astype(float)
        hdu.close()

        flatmodimg=model.flatten()  
        flatresimg=residual.flatten()  

        flatmodimg.sort()
        flatresimg.sort()

        restot=len(flatresimg)

        restop=round(.9*restot)
        resbot=round(.1*restot)

        modimgpatch=flatmodimg
        resimgpatch=flatresimg[resbot:restop]

        modmin = np.min(modimgpatch)
        modmax = np.max(modimgpatch)


        modmin = frac*modmin 
        modmax = fracmax*modmax




        if (modmin > modmax):
            modmin, modmax = modmax, modmin


        resmin = np.min(resimgpatch)
        resmax = np.max(resimgpatch)


        mask=data < 0 
        data[mask] = 1 # avoids problems in log
     
        fig, (ax1, ax2, ax3) = plt.subplots(figsize=(14, 5), nrows = 1, ncols = 3)
        fig.subplots_adjust(left=0.04, right=0.98, bottom=0.02, top=0.98)

        ax1.imshow(con*data + bri, origin ='lower', norm 
                    = colors.LogNorm(vmin=modmin, vmax=modmax), cmap = cmap)

        ax1.set_title('Data')


        for ell in ellipse:
            ax1.add_patch(ell)


        ax2.imshow(con*model + bri, origin='lower', norm 
                    = colors.LogNorm(vmin = modmin, vmax = modmax), cmap = cmap)
        ax2.set_title('GALFIT Model')

        ax3.imshow(residual, origin='lower', vmin = resmin, vmax = resmax, cmap = cmap)
        ax3.set_title('Residual')




        plt.savefig(namepng, dpi = dpival)
    

        #plt.show()



