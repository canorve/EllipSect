
from ellipsect.lib.libs import *

from ellipsect import *



def PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,params,scale):
    """  Produces final best-fitting plot  """

    # subplot for arc sec axis
    plt.close('all')

    # fig, axsec = plt.subplots() #old

    #ULISES begin
    fig, (axsec,axred) = plt.subplots(2, sharex=True, sharey=False)
    gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
    gs.update(hspace=0.07)
    #ULISES end 



    #change the linewidth of the axis
    for axis in ['top','bottom','left','right']:
        axsec.spines[axis].set_linewidth(1.5)



    if params.flagranx[1] == True:
        (xmin,xmax)=params.ranx.split("-")
        xmin=np.float(xmin)
        xmax=np.float(xmax)

    if params.flagrany[1] == True:
        (ymin,ymax)=params.rany.split("-")
        ymin=np.float(ymin)
        ymax=np.float(ymax)


    minrad = np.min(xradq)
    if params.flagranx[1] == False:
        maxrad = np.max(xradq) * params.ranx
    else:
        maxrad = np.max(xradq)

    mincnt = np.min(ysbq)
    maxcnt = np.max(ysbq)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
    yran = mincnt * (maxcnt/mincnt)**np.array([-0.05, +1.05])

    if params.flagrany[1] == False:
        yran1=yran[0]
        yran2=yran[1]

        lyran= yran2 - yran1

        yranmid= yran1 + lyran/2

        lyran=lyran*params.rany

        yran1 = yranmid - lyran/2
        yran2 = yranmid + lyran/2

        yran[0] = yran2 #inverted axis
        yran[1] = yran1



    # ULISES begin
    axsec = plt.subplot(gs[0])
    #axsec.set_xlabel("radius ('')")
    axsec.set_ylabel("Surface Brightness (mag/'')")
    # ULISES end

    axsec.errorbar(xradq, ysbq,yerr=ysberrq,fmt='o-',capsize=2,color='red',markersize=0.7,label="galaxy",linewidth=2)
    axsec.errorbar(xradm, ysbm,yerr=ysberrm,fmt='o-',capsize=2,color='blue',markersize=0.7,label="Model",linewidth=2)
    if params.flagranx[1] == False:
        axsec.set_xlim(xran)
    else:
        axsec.set_xlim(xmin,xmax)


    if params.flagrany[1] == False:
        axsec.set_ylim(yran)
    else:
        axsec.set_ylim(ymax,ymin) #inverted


    if params.flaglogx == True:

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


    # ULISES begin
    axsec.axes.xaxis.set_ticklabels([])
    # Estas dos líneas de abajo las quité mejor porque son cosas muy específicas de mi trabajo así que no tiene caso
    #plt.axvline(x=17, color='k', linestyle='--') 
    #plt.axvline(x=1, color='k', linestyle='dotted')
    # ULISES end 


    axsec.yaxis.set_minor_locator(AutoMinorLocator())
    axsec.yaxis.set_major_locator(AutoLocator())

    if params.flagpix == True:

        axpix = axsec.twiny()

        axpix.set_xlabel("Pixels")
        x1, x2 = axsec.get_xlim()

        axpix.set_xlim(x1/scale, x2/scale)

        axpix.figure.canvas.draw()


        if params.flaglogx == True:
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

    if params.flagrid == True:
        # Customize the major grid
        axsec.grid(which='major', linestyle='-', linewidth='0.7', color='black')
        # Customize the minor grid
        axsec.grid(which='minor', linestyle=':', linewidth='0.5', color='black')

    # ULISES begin
    # Residual plot
    if len(ysbq) < len(ysbm):
        ysbm = ysbm[len(ysbm)-len(ysbq):]
    elif len(ysbq) > len(ysbm):
        ysbq = ysbq[len(ysbq)-len(ysbm):]

    if len(ysbq) < len(ysberrm):
        ysberrm = ysberrm[len(ysberrm)-len(ysbq):]
    elif len(ysbq) > len(ysberrm):
        ysbq = ysbq[len(ysbq)-len(ysberrm):]        
    
    residual = ((ysbq-ysbm)/ysbq)*100 # (data-model)/data in percentage
    err = ysberrm/ysbq*100 # error_model/data in percentage
    axred = plt.subplot(gs[1])
    #axred.scatter(np.log10(xradm),residual, marker='.', color='k')
    if len(xradq) != len(residual):
        axred.errorbar(xradm, residual,yerr=err,fmt='.',capsize=2,color='k')
    else:
        axred.errorbar(xradq, residual,yerr=err,fmt='.',capsize=2,color='k')

    #axsec.errorbar(xradm, ysbm,yerr=ysberrm,fmt='o-',capsize=2,color='blue',markersize=0.7,label="Model",linewidth=2)
    axred.axhline(y=0,ls='dashed', color='k')
    # Estas dos líneas también las quité por lo mismo que dije arriba
    #plt.axvline(x=17, color='k', linestyle='--')
    #plt.axvline(x=1, color='k', linestyle='dotted')
    axred.set_xlabel('Radius (arcsec)')
    axred.set_ylabel('Residual (%)')
    axred.set_ylim(-2,2)
    axred.set_xscale("log")
    # ULISES end



    return xran,yran,axret

#io/plot.py
def PlotSub(xradq,ysbq,nsub,axsec,namec,colorval):
    """
    Produces subcomponent plot

    """

    substr=namec+" "+np.str(nsub+1)

    axsec.plot(xradq, ysbq,'--',color=colorval,linewidth=1.7,markersize=0.7,label=substr)
    #axsec.plot(xradq, ysbq,'--',color=colorval,linewidth=1.5,markersize=0.7,label=substr)


