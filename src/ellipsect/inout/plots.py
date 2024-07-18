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

from ellipsect.sectors.num import Interpol
from ellipsect.inout.prt import  PrintFilesGax
from ellipsect.inout.prt import PrintFilesComps

from ellipsect.inout.galfit  import conver2Sersic 
from ellipsect.inout.galfit  import GetReff 

def PlotSB(xradq, ysbq, ysberrq, xradm, ysbm, ysberrm, ellconf, scale):
    """  Produces final best-fitting plot  """

    # subplot for arc sec axis
    plt.close('all')



    #ULISES begin
    #fig, (axsec, axred) = plt.subplots(2, sharex=True, sharey=False)
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

    axsec.set_ylabel(r"Surface Brightness $(mag\; arcsec^{-2})$")
    # ULISES end

    axsec.errorbar(xradq, ysbq, yerr=ysberrq,fmt='o-',capsize=2,color='red',markersize=0.7,label="galaxy",linewidth=2)
    #axsec.plot(r, mgegal.sb[angal], 'C3-',linewidth=2)


    if ellconf.flagalax == False:
        axsec.errorbar(xradm, ysbm, yerr=ysberrm,fmt='o-',capsize=2,color='blue',markersize=0.7,label="Model",linewidth=2)


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
        axsec.axvline(x=xpos,  linestyle='--', color='k', label='FWHM', linewidth=2)
    # end 

    #effective radius 
    if ellconf.flagrep: 
        xre = ellconf.rep*scale
        axsec.axvline(x=xre,  linestyle='--', color='r', label='Re', linewidth=1.5)
 
    #90% of total light radius 
    if ellconf.flagr90p: 
        xr90 = ellconf.r90p*scale
        axsec.axvline(x = xr90,  linestyle='--', color='b', label='R90', linewidth=1.5)
 
    #95%  of total light radius 
    if ellconf.flagr95p: 
        xr95 = ellconf.r95p*scale
        axsec.axvline(x = xr95,  linestyle='--', color='c', label='R95', linewidth=1.5)
 

    # ULISES begin
    axsec.axes.xaxis.set_ticklabels([])
    # ULISES end 


    axsec.yaxis.set_minor_locator(AutoMinorLocator())
    axsec.yaxis.set_major_locator(AutoLocator())

    if ellconf.title:
        plt.title(ellconf.namefile, fontsize=10)

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
    axred.set_xlabel('Radius $(arcsec)$')
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



    return xran, yran, axret

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
              cmap = ellconf.cmap, ellipse = ell, plate = galhead.scale)


    if ellconf.dplot:
        plt.pause(1.5)
 
    plt.close()

    #}
    #####################


def Comp2Ellip(galhead, galcomps, N, lw=1):
    ''' converts galfit component parameter into an Ellipse object''' 


    ellipses = [] 

    #color value
    values = np.arange(N)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


    idxcol=0
    for idx, item in enumerate(galcomps.N):

        if galcomps.Active[idx] == True:
            # correcting coordinates
            xc = galcomps.PosX[idx] - galhead.xmin + 1
            yc = galcomps.PosY[idx] - galhead.ymin + 1


            pa = galcomps.PosAng[idx] + 90

            w = galcomps.Rad[idx]
            h = galcomps.Rad[idx]*galcomps.AxRat[idx]

            colorVal = scalarMap.to_rgba(values[idxcol])


            ell=Ellipse((xc, yc), width = w, height = h, angle = pa,
                         edgecolor = colorVal,
                         facecolor = 'none',
                         linewidth = lw)

            ellipses.append(ell)
            idxcol+=1

    return ellipses






class ShowCube:

    def __init__(self, cubeimg: str, namepng="cubeout.png", dpival=100, 
                bri = 0, con = 1, cmap='viridis', ellipse=[], plate = 1):
        """
        This routine shows the GALFIT output cube image: galaxy, model and residual    
        """

        
        hdu = fits.open(cubeimg)
        data = (hdu[1].data.copy()).astype(float)
        model = (hdu[2].data.copy()).astype(float)
        residual = (hdu[3].data.copy()).astype(float)
        hdu.close()

        flatmodimg = model.flatten()  
        flatresimg = residual.flatten()  

        flatmodimg.sort()
        flatresimg.sort()

        restot = len(flatresimg)

        restop = round(.9*restot)
        resbot = round(.1*restot)

        modtot = len(flatmodimg)

        modtop = round(.9*modtot)
        modbot = round(.1*modtot)



        modimgpatch = flatmodimg#[modbot:modtop]
        resimgpatch = flatresimg[resbot:restop]

        resmin = np.min(resimgpatch)
        resmax = np.max(resimgpatch)


        modmin = np.min(modimgpatch)
        modmax = np.max(modimgpatch)


        data = data.clip(modmax/1e4,modmax)
        model = model.clip(modmax/1e4)

        modmin = modmax/1e4


        middle = (modmax - modmin)/2


        #brightness auto-adjust according to the contrast value 

        Autobri = middle*(con -1) + modmin*(1-con) 


        #user can re-adjust brightness in addition to Autobri
        newdata = con*(data - middle) + middle + Autobri + bri*(modmax-middle)
        newmodel = con*(model - middle) + middle + Autobri + bri*(modmax-middle)



        mask=data < 0 
        data[mask] = 1 # avoids problems in log
     
        fig, (ax1, ax2, ax3) = plt.subplots(figsize=(14, 5), nrows = 1, ncols = 3)

        im = ax1.imshow(newdata, origin ='lower', interpolation='nearest', norm 
                    = colors.LogNorm(vmin=modmin, vmax=modmax), cmap = cmap)


        ax1.set_title('Data')

        y,x = data.shape

        xt = .02*x
        yt = .02*y

        lxline = round(.1*x)
        lyline = round(.1*y)
        x1 = [xt, xt+lxline]
        y1 = [lyline, lyline]

        arcsec = lxline*plate*U.arcsec 


        if arcsec.value >= 60:
            lxlinearc = arcsec.to("arcmin").value
            s = "{}\'".format(round(lxlinearc))
        else:
            lxlinearc = arcsec.value
            s = "{}\'\'".format(round(lxlinearc))
 
        ax1.plot(x1, y1, color="white", linewidth=3)

        ax1.text(xt+round(lxline/5),lyline+yt,s,color='white',fontsize=14)
    
        #ax1.set_xlabel(r'$\circ$')
        #ax2.set_ylabel('23\"')

        for ell in ellipse:
            ax1.add_patch(ell)


        ax2.imshow(newmodel, origin='lower', interpolation='nearest', norm 
                    = colors.LogNorm(vmin = modmin, vmax = modmax), cmap = cmap)


        ax2.set_title('GALFIT Model')

        ax3.imshow(residual, origin='lower', vmin = resmin, vmax = resmax, cmap = cmap)
        ax3.set_title('Residual')


        fig.subplots_adjust(left=0.04, right=0.98, bottom=0.02, top=0.98)
        plt.savefig(namepng, dpi = dpival)
    

        #plt.show()


def PlotMul(ellconf, galhead, galcomps, mgegal, mgemod, mgecom):


    
    comps = conver2Sersic(galcomps) #to compute Re


    fignum = 1

    minrad = np.min(mgegal.rad)
    maxrad = np.max(mgegal.rad)

    minsb = np.min(mgegal.sb)
    maxsb = np.max(mgegal.sb)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
    yran = minsb * (maxsb/minsb)**np.array([+1.05,-0.05]) #inverted axis


    sectors = np.unique(mgegal.angle)
    n = sectors.size
    dn = int(round(n/6.))
    nrows = (n-1)//dn + 1 # integer division




    if ellconf.flagranx == True:
        (xmin,xmax) = ellconf.ranx[0], ellconf.ranx[1]

    if ellconf.flagrany == True:
        (ymin,ymax) = ellconf.rany[0], ellconf.rany[1]



    plt.clf()


    if ellconf.title:
        plt.title(ellconf.namefile, fontsize=10)



    fig, axsec = plt.subplots(nrows, 2, sharex=True, sharey='col', num=fignum)
    fig.subplots_adjust(hspace=0.01)






    if ellconf.flagpix:
        axpix = axsec[0,0].twiny()
        axpix2 = axsec[0,1].twiny()

    fig.text(0.04, 0.5, r"Surface Brightness $(mag\; arcsec^{-2})$", va='center', rotation='vertical')
    fig.text(0.96, 0.5, 'error ', va='center', rotation='vertical')

    axsec[-1, 0].set_xlabel("radius $(arcsec)$")
    axsec[-1, 1].set_xlabel("radius $(arcsec)$")

    if ellconf.flaglogx == True:
        axsec[-1, 0].xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        if ellconf.flagpix:
            axpix.set_xscale("log")
            axpix2.set_xscale("log")

    else:
        axsec[-1, 0].xaxis.set_minor_locator(AutoMinorLocator())

    axsec[-1, 0].tick_params(which='both', width=2)
    axsec[-1, 0].tick_params(which='major', length=7)
    axsec[-1, 0].tick_params(which='minor', length=4, color='r')

    if ellconf.flaglogx == True:
        axsec[-1, 0].xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
    else:
        axsec[-1, 0].xaxis.set_minor_locator(AutoMinorLocator())
    axsec[-1, 1].tick_params(which='both', width=2)
    axsec[-1, 1].tick_params(which='major', length=7)
    axsec[-1, 1].tick_params(which='minor', length=4, color='r')

    angsec = 90 - ellconf.parg

    #    row = 7 # old values
    #row = nrows -1
    row = 0  # major axis in first row

    for j in range(0, n, dn):
        angal = np.nonzero(mgegal.angle == sectors[j])[0]

        angal = angal[np.argsort(mgegal.rad[angal])]
        r = mgegal.rad[angal]

        angmod = np.nonzero(mgemod.angle == sectors[j])[0]
        angmod = angmod[np.argsort(mgemod.rad[angmod])]

        #if (len(mgemodrad) < len(mgerad)):
        #    r2 = mgemodrad[angmod]
        #else:
        #    angmod=w
        #    r2 = mgemodrad[angmod]

        r2 = mgemod.rad[angmod]

        #angsec=90-ellconf.parg
        txtang = sectors[j]
        txtangsky = sectors[j] + ellconf.parg #angle measured from sky north. Same as GALFIT

        if txtangsky > 90:
            txtangsky = txtangsky - 180 


        txt = r"$%.f^\circ$" % txtang
        txtsky = r"$%.f^\circ$" % txtangsky

        txtminor= "minor axis"
        txtmajor= "major axis"

        if ellconf.flagranx == True:
            axsec[row, 0].set_xlim(xmin, xmax)
        else:
            axsec[row, 0].set_xlim(xran)

        if ellconf.flagrany == True:
            axsec[row, 0].set_ylim(ymax, ymin) #inverted
        else:
            axsec[row, 0].set_ylim(yran)


        #begin psf fwhm 
        if ellconf.flagfwhm: 
            xpos = ellconf.fwhm*galhead.scale
            axsec[row, 0].axvline(x=xpos, linestyle='--', color='k', linewidth=2)
        # end 


        #effective radius 
        if ellconf.flagrep: 
            eff = 0.5
            xre, tempmag = GetReff().GetReSer(galhead, comps, eff, txtangsky)

            xre = xre*galhead.scale
            axsec[row, 0].axvline(x=xre,  linestyle='--', color='r', label='Re', linewidth=1.5)
     
        #90% radius 
        if ellconf.flagr90p: 
            eff = 0.9
            xr90, tempmag = GetReff().GetReSer(galhead, comps, eff, txtangsky)

            xr90 = xr90*galhead.scale
            axsec[row, 0].axvline(x = xr90,  linestyle='--', color='b', label='R90', linewidth=1.5)
 
        #95% radius 
        if ellconf.flagr95p: 
            eff = 0.95
            xr95, tempmag = GetReff().GetReSer(galhead, comps, eff, txtangsky)

            xr95 = xr95*galhead.scale
            axsec[row, 0].axvline(x = xr95,  linestyle='--', color='c', label='R95', linewidth=1.5)
 



        if ellconf.flaglogx == False:

            #axsec[row, 0].plot(r, mgesb[angal], 'C3o') 
            #change lines instead of dots
            axsec[row, 0].plot(r, mgegal.sb[angal], 'C3-',linewidth=2)

            if ellconf.flagalax == False:
                axsec[row, 0].plot(r2, mgemod.sb[angmod], 'C0-', linewidth=1.5)

        else:

            #axsec[row, 0].semilogx(r, mgesb[angal], 'C3o')
            #change lines instead of dots
            axsec[row, 0].semilogx(r, mgegal.sb[angal], 'C3-', linewidth=2)

            if ellconf.flagalax == False:
                axsec[row, 0].semilogx(r2, mgemod.sb[angmod], 'C0-', linewidth=1.5)

        if ellconf.flagsbout == True: 

            rtxtang=np.int32(np.round(txtang)) 

            PrintFilesGax(ellconf,galhead,rtxtang,r,mgegal.sb,angal,r2,mgemod.sb,angmod)


        if ellconf.flagrid == True:
            # Customize the major grid
            axsec[row,0].grid(which='major', linestyle='-', linewidth='0.7', color='black')
            # Customize the minor grid
            axsec[row,0].grid(which='minor', linestyle=':', linewidth='0.5', color='black')

            #  axsec[row,0].grid(True)

        maskgal = galcomps.Active == True 
        if ellconf.flagcomp == True:
            ii=0
                #color value

            values = range(len(galcomps.N[maskgal]))
            jet = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

            while(ii<len(galcomps.N[maskgal])):

                #angtemp = np.nonzero(mgeanglesub[ii] == sectors[j])[0]
        
                ######################## Patch for angle :############  
                alpha = sectors[j]
                angsec2 = 90 - galcomps.PosAng[maskgal][ii]

                if angsec < 0:
                    angsec = 360 + angsec
                if angsec2 < 0:
                    angsec2 = 360 + angsec2

                alpha2 = alpha  + angsec - angsec2 

                if alpha2 < 0:
                    alpha2 = 360 + alpha2

                if alpha2 > 90 and alpha2 <=270:
                    alpha2 = np.abs(180-alpha2)

                if alpha2 > 270 and alpha2 <=360:
                    alpha2 = np.abs(360-alpha2)


                # search for the nearest angle for subcomponent:
                jj=(np.abs(mgecom.sector[ii] - alpha2)).argmin()  

                diffangle =  mgecom.sector[ii][jj] - alpha2

                # alpha: angle from major axis of galaxy
                # angsec: position angle of the galaxy
                # theta2: position angle of the component
                # alpha2: angle from major axis of component
                # sectors: angle obtained from sectors_photometry
                # it is expected that alpha2 and sectors are the closest possible.

                ###############################################
                
                #angtemp = np.nonzero(mgeanglesub[ii] == sectorsub[ii][j])[0]
                angtemp = np.nonzero(mgecom.angle[ii] == mgecom.sector[ii][jj])[0]
                angtemp = angtemp[np.argsort(mgecom.rad[ii][angtemp])]

      
                rtemp = mgecom.rad[ii][angtemp]

                colorval = scalarMap.to_rgba(values[ii])
                if ellconf.flaglogx == False:
                #    axsec[row, 0].plot(rtemp, mgesbsub[ii][angtemp],'--',color='skyblue', linewidth=2)
                    axsec[row, 0].plot(rtemp, mgecom.sb[ii][angtemp],'--',color=colorval, linewidth=1.5)
                else:
                    axsec[row, 0].semilogx(rtemp, mgecom.sb[ii][angtemp], '--',color=colorval, linewidth=1.5)

                if ellconf.flagsbout == True:
                    ncomp=ii+1
                    ncomp=str(ncomp)

                    PrintFilesComps(ellconf,galhead,galcomps,rtxtang,ncomp,
                                    diffangle,rtemp,mgecom.sb,ii,angtemp)

                ii+=1

        #axsec[row, 0].text(0.98, 0.95, txt, ha='right', va='top', transform=axsec[row, 0].transAxes)
        axsec[row, 0].text(0.98, 0.95, txtsky, color='red',ha='right', va='top', transform=axsec[row, 0].transAxes)
        axsec[row, 0].text(0, 0, txt, ha='left', va='bottom', color='grey', transform=axsec[row, 0].transAxes)

        if (len(mgemod.rad) > len(mgegal.rad)):

            mgemodsbnew, smooth_flag = Interpol(r2,mgemod.sb[angmod],r)
            sberr = 1 - mgemodsbnew/mgegal.sb[angal]
            axsec[row, 1].plot(r, sberr*100, 'C0o')
            if(smooth_flag):
                print("smoothing interpolation was used for angle: ",np.int32(np.round(txtang)))


        else:
            mgesbnew, smooth_flag = Interpol(r,mgegal.sb[angal],r2)
            sberr=1-mgemod.sb[angmod]/mgesbnew
            axsec[row, 1].plot(r2, sberr*100, 'C0o')

            if(smooth_flag):
                print("smoothing interpolation was used for angle: ",np.int32(np.round(txtang)))

        axsec[row, 1].axhline(linestyle='--', color='C1', linewidth=2)
        axsec[row, 1].yaxis.tick_right()
        axsec[row, 1].yaxis.set_label_position("right")
        axsec[row, 1].set_ylim([-19.5, 20])
        # axsec[row, 1].set_ylim([-20, 20])
        #axsec[row, 1].text(0.98, 0.95, txt, ha='right', va='top', transform=axsec[row, 1].transAxes)
        axsec[row, 1].text(0.98, 0.95, txtsky, fontweight='bold', 
                            color='red',ha='right', va='top', 
                            transform = axsec[row, 1].transAxes)

        axsec[row, 1].text(0, 0, txt, ha='left', va='bottom',
                            color='grey', transform=axsec[row, 1].transAxes)

        if (txtang == 0):
            axsec[row, 1].text(0.98, 0.10, txtmajor, fontweight='bold',
                                fontsize=8.5, ha='right', va='bottom', 
                                transform = axsec[row, 1].transAxes)

        if (txtang == 90):
            axsec[row, 1].text(0.98, 0.10, txtminor, fontweight='bold',
                                fontsize=8.5, ha='right', va='bottom', 
                                transform=axsec[row, 1].transAxes)


        if ellconf.flagranx == True:
            axsec[row, 1].set_xlim(xmin,xmax)
        else:
            axsec[row, 1].set_xlim(xran)

        axsec[row, 0].yaxis.set_minor_locator(AutoMinorLocator())
        axsec[row, 0].tick_params(which='both', width=2)
        axsec[row, 0].tick_params(which='major', length=7)
        axsec[row, 0].tick_params(which='minor', length=4, color='r')

        axsec[row, 1].yaxis.set_minor_locator(AutoMinorLocator())
        axsec[row, 1].tick_params(which='both', width=2)
        axsec[row, 1].tick_params(which='major', length=7)
        axsec[row, 1].tick_params(which='minor', length=4, color='r')


        #change the linewidth of the axis
        for axis in ['top','bottom','left','right']:
            axsec[row,0].spines[axis].set_linewidth(1.5)
            axsec[row,1].spines[axis].set_linewidth(1.5)


        # row -= 1
        row += 1


    if ellconf.flagpix == True:
        axpix.set_xlabel("(pixels)")
     
        #x1, x2 = axsec[7,0].get_xlim() ## buggy for some data have to change it for code below:
        
        if ellconf.flagranx == True:
            x1=xmin
            x2=xmax
        else:
            x1= xran[0]
            x2= xran[1]
 
        axpix.set_xlim(x1/galhead.scale, x2/galhead.scale)
        axpix.figure.canvas.draw()

        axpix2.set_xlabel("(pixels)")
        axpix2.set_xlim(x1/galhead.scale, x2/galhead.scale)
        axpix2.figure.canvas.draw()

        ##
        if ellconf.flaglogx == True:
            axpix.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        else:
            axpix.xaxis.set_minor_locator(AutoMinorLocator())

        axpix.tick_params(which='both', width=2)
        axpix.tick_params(which='major', length=7)
        axpix.tick_params(which='minor', length=4, color='r')

        if ellconf.flaglogx == True:
            axpix2.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        else:
            axpix2.xaxis.set_minor_locator(AutoMinorLocator())

        axpix2.tick_params(which='both', width=2)
        axpix2.tick_params(which='major', length=7)
        axpix2.tick_params(which='minor', length=4, color='r')







