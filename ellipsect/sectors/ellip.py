

from ellipsect.lib.libs import *


#next function removed from this file and moved to comp
from ellipsect.sectors.comp import FindSB
from ellipsect.inout.plots import PlotSB

from ellipsect.sectors.num import Interpol
from ellipsect.sectors.num import Re90

from ellipsect.sectors.comp import SubComp


from ellipsect.inout.prt import PrintEllFilesGax
from ellipsect.inout.prt import PrintFilesGax



def EllipSectors(params, galpar, galcomps, sectgalax, sectmodel, sectcomps,n_sectors=19, minlevel=0):


    xradm = []
    ysbm = []
    ysberrm = []


    # galax 

    stidxg = np.argsort(sectgalax.radius)

    mgerad=sectgalax.radius[stidxg]
    mgecount=sectgalax.counts[stidxg]
    mgeangle=sectgalax.angle[stidxg]
    mgeanrad=np.deg2rad(mgeangle)

    ab=galpar.q

    aellabg= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)

    #changing to arc sec
    aellarcg=aellabg*galpar.scale


    ####
    # formula according to cappellary mge manual:
    mgesbg= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
    ####

    stidxq = np.argsort(aellarcg)


    xarcg = aellarcg[stidxq]
    ymgeg = mgesbg[stidxq]

    ymgec = mgecount[stidxq] + galpar.skylevel

    #############  Function to order SB along X-axis for galaxy

    #xradq, ysbq, ysberrq    = FindSB(xarcg, ymgeg, n_sectors)

    xradq, ysbq, ysberrq, ysbc, ysbcerr  = FindSBCounts(xarcg, ymgeg, ymgec, n_sectors)


    #######################################
    #######################################
    #######################################

    # computing sky as a reference.
    if params.flagradsky:

        # computing sky with the gradient method
        print("Computing sky as a reference. This will not be used for output computations.")

        ImageFile = galpar.inputimage
        MaskFile = galpar.maskimage

        xx = galpar.inxc
        yy = galpar.inyc

        thetadeg = galpar.ang
        #e = 1 - galpar.q
        q = galpar.q

        width = params.skywidth


        ###
        Rinit = 1

        if not(params.flagskyRad):
            #rad90= galpar.rad * (1.53 + 0.73 * galpar.serind+ 0.07 * galpar.serind**2) 
            rad90= Re90(galpar.rad,galpar.serind) 
            Rinit = 1*rad90 # 1 times the R 90% of light radius

            if (Rinit < 50): # Rinit can not be less than the default value 
                Rinit = 50 
            params.skyRad = Rinit # save value for output

        else:
            Rinit = params.skyRad
            


        print("computing sky with the gradient method")

        line="using Rinit = {:.2f} width = {}".format(Rinit,width)
        print(line)

        line="using thetadeg = {:.2f} q = {}".format(thetadeg,galpar.q)
        print(line)
        
        line="using xx = {} yy  = {}".format(xx,yy)
        print(line)

        mean,std, median,rad = SkyCal().GetEllipSky(ImageFile,MaskFile,xx,yy,thetadeg,q,Rinit,width,params.namering)

        line="Total sky:  mean = {:.2f}; std={:.2f}; median = {} ".format(mean,std,median)
        print(line)

        #saving for output
        galpar.gradskymean = mean  
        galpar.gradskystd = std
        galpar.gradskymed = median





    if params.flagrandboxsky:

        # computing sky  using random boxes across the image
        print("Computing sky as a reference. This will not be used for output computations.")

        ImageFile = galpar.inputimage
        MaskFile = galpar.maskimage

        xx = galpar.xc
        yy = galpar.yc

        thetadeg = galpar.ang
        #e = 1 - galpar.q
        q = galpar.q

        ###

        box = params.skybox
        num = params.skynum


        ###
        Rinit = 1

        if not(params.flagskyRad):
            #rad90=galpar.rad * (1.53 + 0.73 * galpar.serind + 0.07 * galpar.serind**2) 
            rad90= Re90(galpar.rad,galpar.serind) 
            Rinit = 3*rad90 # 1 times the R 90% of light radius

            if (Rinit < 100): # Rinit can not be less than the default value 
                Rinit = 100 
            params.skyRad = Rinit # save value for output
        else:
            Rinit = params.skyRad


        print("computing sky with the random box method")

        line="using Rad = {:.2f}, box size= {}, number of boxes = {}".format(Rinit,box,num)
        print(line)

        ##
        if params.flagskyRadmax:
            Rmax = params.skyRadmax
            mean,std, median = SkyCal().RandBox(ImageFile,MaskFile,xx,yy,thetadeg,q,Rinit,box,num,Rmax)
        else:
            Rmax = 0
            mean,std, median = SkyCal().RandBox(ImageFile,MaskFile,xx,yy,thetadeg,q,Rinit,box,num,Rmax)
        #

        line="Total sky:  mean = {:.2f}; std = {:.2f}; median = {:.2f}".format(mean,std,median)
        print(line)

        #saving for output
        galpar.randskymean = mean  
        galpar.randskystd = std
        galpar.randskymed = median




    #######################################
    #######################################
    #######################################

    # model
    stidxm = np.argsort(sectmodel.radius)

    mgerad=sectmodel.radius[stidxm]
    mgecount=sectmodel.counts[stidxm]
    mgeangle=sectmodel.angle[stidxm]
    mgeanrad=np.deg2rad(mgeangle)

    ab=galpar.q

    aellabm= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)


    aellarcm=aellabm*galpar.scale

    # formula according to cappellary mge manual
    mgesbm= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
    ##

    stidxq = np.argsort(aellarcm)

    xarcm = aellarcm[stidxq]
    ymgem = mgesbm[stidxq]




    ######  Function to order SB along X-axis for model

    xradm, ysbm, ysberrm    = FindSB(xarcm, ymgem, n_sectors)

    ################ Plotting

    limx,limy,axsec=PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,params,galpar.scale)


    ### surface brightness output file

    if params.flagsbout == True: 

        #print to file    
        PrintEllFilesGax(params,galpar,xradq,ysbq,ysberrq,xradm,ysbm,ysberrm)


    #### Creating Subcomponents images with Galfit


    if params.flagcomp:


        xradq,ysbq,n=SubComp(params, galpar, galcomps, sectcomps, axsec, n_sectors=n_sectors)



    axsec.legend(loc=1)

    return limx,limy

#sectors/ellipsectors.py
def MulEllipSectors(params, galpar, galcomps, sectgalax, sectmodel, sectcomps):

  
    fignum=1


    eps=1-galpar.q

    if params.flagranx[1] == True:
        (xmin,xmax)=params.ranx.split("-")
        xmin=np.float(xmin)
        xmax=np.float(xmax)

    if params.flagrany[1] == True:
        (ymin,ymax)=params.rany.split("-")
        ymin=np.float(ymin)
        ymax=np.float(ymax)


    yctemp=galpar.xc
    xctemp=galpar.yc


    angsec=90-galpar.ang

    ######################

    sg = sectgalax

    sm = sectmodel
    ###################################################

    stidx = np.argsort(sg.radius)

    #   galaxy
    mgerad=sg.radius[stidx]

    mgecount=sg.counts[stidx]
    mgeangle=sg.angle[stidx]
    mgeanrad=np.deg2rad(mgeangle)


    # model

    stidx = np.argsort(sm.radius)

    mgemodrad=sm.radius[stidx]

    mgemodcount=sm.counts[stidx]
    mgemodangle=sm.angle[stidx]
    mgemodanrad=np.deg2rad(mgemodangle)


    # converting to pixels

    mgerad=mgerad*galpar.scale
    mgemodrad=mgemodrad*galpar.scale


    # formula according to cappellary mge manual
    # galaxy:
    mgesb= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
    # Model:
    mgemodsb= galpar.mgzpt - 2.5*np.log10(mgemodcount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1


    if params.flagcomp:

        angtemp=[]
        mgesbsub=[]
        mgeradsub=[]
        mgeanglesub=[]
        sectorsub=[]

        ###############################
        
        ab=galpar.q
        ni=0
        while(ni<len(galcomps.N)):

            subcmp = sectcomps[ni]


            subidx = np.argsort(subcmp.radius)

            temprad=subcmp.radius[subidx]

            #converting to arcsec
            temprad=temprad*galpar.scale

            mgecountsub=subcmp.counts[subidx]

            tempangle=subcmp.angle[subidx]
            mgeanradsub=np.deg2rad(tempangle)


            # formula according to cappellary mge manual
            tempmge= galpar.mgzpt - 2.5*np.log10(mgecountsub/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
            

            tempsectorsub = np.unique(tempangle)

            sectorsub.append(tempsectorsub)
            mgesbsub.append(tempmge)
            mgeradsub.append(temprad)
            mgeanglesub.append(tempangle)

            ni+=1


    minrad = np.min(mgerad)

    if params.flagranx[1] == False:
        maxrad = np.max(mgerad) * params.ranx
    else:
        maxrad = np.max(mgerad)

    minsb = np.min(mgesb)
    maxsb = np.max(mgesb)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
    yran = minsb * (maxsb/minsb)**np.array([+1.05,-0.05])


    if params.flagrany[1] == False:

        yran1=yran[0]
        yran2=yran[1]

        lyran= yran2 - yran1

        yranmid= yran1 + lyran/2

        lyran=lyran*params.rany

        yran1 = yranmid - lyran/2
        yran2 = yranmid + lyran/2

        yran[0] = yran1
        yran[1] = yran2


    sectors = np.unique(mgeangle)
    n = sectors.size
    dn = int(round(n/6.))
    nrows = (n-1)//dn + 1 # integer division

    plt.clf()

    fig, axsec = plt.subplots(nrows, 2, sharex=True, sharey='col', num=fignum)
    fig.subplots_adjust(hspace=0.01)


    if params.flagpix:
        axpix = axsec[0,0].twiny()
        axpix2 = axsec[0,1].twiny()

    fig.text(0.04, 0.5, "Surface Brightness (mag/'')", va='center', rotation='vertical')
    fig.text(0.96, 0.5, 'error (%)', va='center', rotation='vertical')

    axsec[-1, 0].set_xlabel("radius ('')")
    axsec[-1, 1].set_xlabel("radius ('')")

    if params.flaglogx == True:
        axsec[-1, 0].xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        if params.flagpix:
            axpix.set_xscale("log")
            axpix2.set_xscale("log")

    else:
        axsec[-1, 0].xaxis.set_minor_locator(AutoMinorLocator())

    axsec[-1, 0].tick_params(which='both', width=2)
    axsec[-1, 0].tick_params(which='major', length=7)
    axsec[-1, 0].tick_params(which='minor', length=4, color='r')

    if params.flaglogx == True:
        axsec[-1, 0].xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
    else:
        axsec[-1, 0].xaxis.set_minor_locator(AutoMinorLocator())
    axsec[-1, 1].tick_params(which='both', width=2)
    axsec[-1, 1].tick_params(which='major', length=7)
    axsec[-1, 1].tick_params(which='minor', length=4, color='r')


    #    row = 7 # old values
    #row = nrows -1
    row = 0  # major axis in first row

    for j in range(0, n, dn):
        angal = np.nonzero(mgeangle == sectors[j])[0]

        angal = angal[np.argsort(mgerad[angal])]
        r = mgerad[angal]

        angmod = np.nonzero(mgemodangle == sectors[j])[0]
        angmod = angmod[np.argsort(mgemodrad[angmod])]

        #if (len(mgemodrad) < len(mgerad)):
        #    r2 = mgemodrad[angmod]
        #else:
        #    angmod=w
        #    r2 = mgemodrad[angmod]

        r2 = mgemodrad[angmod]

        txtang= sectors[j]
        txt = r"$%.f^\circ$" % txtang

        if params.flagranx[1] == False:
            axsec[row, 0].set_xlim(xran)
        else:
            axsec[row, 0].set_xlim(xmin,xmax)

        if params.flagrany[1] == False:
            axsec[row, 0].set_ylim(yran)
        else:
            axsec[row, 0].set_ylim(ymax,ymin) #inverted

        if params.flaglogx == False:

            #axsec[row, 0].plot(r, mgesb[angal], 'C3o') 
            #change lines instead of dots
            axsec[row, 0].plot(r, mgesb[angal], 'C3-',linewidth=2)

            axsec[row, 0].plot(r2, mgemodsb[angmod], 'C0-', linewidth=1.5)

        else:

            #axsec[row, 0].semilogx(r, mgesb[angal], 'C3o')
            #change lines instead of dots
            axsec[row, 0].semilogx(r, mgesb[angal], 'C3-', linewidth=2)

            axsec[row, 0].semilogx(r2, mgemodsb[angmod], 'C0-', linewidth=1.5)

        if params.flagsbout == True: 

            rtxtang=np.int(np.round(txtang)) 

            PrintFilesGax(params,galpar,rtxtang,r,mgesb,angal,r2,mgemodsb,angmod)


        if params.flagrid == True:
            # Customize the major grid
            axsec[row,0].grid(which='major', linestyle='-', linewidth='0.7', color='black')
            # Customize the minor grid
            axsec[row,0].grid(which='minor', linestyle=':', linewidth='0.5', color='black')

            #  axsec[row,0].grid(True)

        if params.flagcomp == True:
            ii=0
                #color value
            values = range(len(galcomps.N))
            jet = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

            while(ii<len(galcomps.N)):

                #angtemp = np.nonzero(mgeanglesub[ii] == sectors[j])[0]
        
                ######################## Patch for angle :############  
                alpha = sectors[j]
                angsec2= 90-galcomps.PosAng[ii]
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
                jj=(np.abs(sectorsub[ii]-alpha2)).argmin()  

                diffangle =  sectorsub[ii][jj] - alpha2
                # uncomment below to check angles in multiplot:
                #print("Check C{}: Axrat: {:.3f}, alpha: {:.3f} angsec: {:.3f} ; theta2: {:.3f} sector {:.3f}; alpha2 {:.3f} ".format(ii,galcomps.AxRat[ii],alpha,angsec,90-galcomps.PosAng[ii],sectorsub[ii][jj],alpha2))
                # alpha: angle from major axis of galaxy
                # angsec: position angle of the galaxy
                # theta2: position angle of the component
                # alpha2: angle from major axis of component
                # sectors: angle obtained from sectors_photometry
                # it is expected that alpha2 and sectors are the closest possible.

                ###############################################
                
                #angtemp = np.nonzero(mgeanglesub[ii] == sectorsub[ii][j])[0]
                angtemp = np.nonzero(mgeanglesub[ii] == sectorsub[ii][jj])[0]
                angtemp = angtemp[np.argsort(mgeradsub[ii][angtemp])]

      
                rtemp = mgeradsub[ii][angtemp]

                colorval = scalarMap.to_rgba(values[ii])
                if params.flaglogx == False:
                #    axsec[row, 0].plot(rtemp, mgesbsub[ii][angtemp],'--',color='skyblue', linewidth=2)
                    axsec[row, 0].plot(rtemp, mgesbsub[ii][angtemp],'--',color=colorval, linewidth=1.5)
                else:
                    axsec[row, 0].semilogx(rtemp, mgesbsub[ii][angtemp], '--',color=colorval, linewidth=1.5)

                #introduce output 
                if params.flagsbout == True:
                    ncomp=ii+1
                    ncomp=str(ncomp)

                    PrintFilesComps(params,galpar,galcomps,rtxtang,ncomp,diffangle,rtemp,mgesbsub,ii,angtemp)

                ii+=1

        axsec[row, 0].text(0.98, 0.95, txt, ha='right', va='top', transform=axsec[row, 0].transAxes)


        if (len(mgemodrad) > len(mgerad)):

            mgemodsbnew = Interpol(r2,mgemodsb[angmod],r)
            sberr=1-mgemodsbnew/mgesb[angal]
            axsec[row, 1].plot(r, sberr*100, 'C0o')

        else:
            mgesbnew = Interpol(r,mgesb[angal],r2)
            sberr=1-mgemodsb[angmod]/mgesbnew
            axsec[row, 1].plot(r2, sberr*100, 'C0o')


        axsec[row, 1].axhline(linestyle='--', color='C1', linewidth=2)
        axsec[row, 1].yaxis.tick_right()
        axsec[row, 1].yaxis.set_label_position("right")
        axsec[row, 1].set_ylim([-19.5, 20])
        # axsec[row, 1].set_ylim([-20, 20])
        axsec[row, 1].text(0.98, 0.95, txt, ha='right', va='top', transform=axsec[row, 1].transAxes)


        if params.flagranx[1] == False:
            axsec[row, 1].set_xlim(xran)
        else:
            axsec[row, 1].set_xlim(xmin,xmax)


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


    if params.flagpix == True:
        axpix.set_xlabel("(pixels)")
     
        #x1, x2 = axsec[7,0].get_xlim() ## buggy for some data have to change it for code below:
        
        if params.flagranx[1] == False:
            x1= xran[0]
            x2= xran[1]
        else:
            x1=xmin
            x2=xmax

        axpix.set_xlim(x1/galpar.scale, x2/galpar.scale)
        axpix.figure.canvas.draw()

        axpix2.set_xlabel("(pixels)")
        axpix2.set_xlim(x1/galpar.scale, x2/galpar.scale)
        axpix2.figure.canvas.draw()

        ##
        if params.flaglogx == True:
            axpix.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        else:
            axpix.xaxis.set_minor_locator(AutoMinorLocator())

        axpix.tick_params(which='both', width=2)
        axpix.tick_params(which='major', length=7)
        axpix.tick_params(which='minor', length=4, color='r')

        if params.flaglogx == True:
            axpix2.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        else:
            axpix2.xaxis.set_minor_locator(AutoMinorLocator())
        axpix2.tick_params(which='both', width=2)
        axpix2.tick_params(which='major', length=7)
        axpix2.tick_params(which='minor', length=4, color='r')


#sectors/ellipsectors.py
def FindSBCounts(xarcq, ymgeq, ymgec, numsectors):
    # the xarcq array must be ordered
    # use mag instead of counts

    xradq=[]
    ysbq=[]
    ysberrq=[]
    ysbc=[]
    ysbcerr=[]

    xradq=np.array(xradq)
    ysbq=np.array(ysbq)
    ysberrq=np.array(ysberrq)
    ysbc=np.array(ysbc)
    ysbcerr=np.array(ysbcerr)


    numsave=0
    tot=xarcq.size
    count=0
    for i in range(tot,0,-1):

        lima=i-numsectors
        limb=i

        if xarcq[lima:limb].size == 0:
            break
        else:
            valstd=np.std(xarcq[lima:limb])
            if valstd < 0.1:
                numsave=count
                break
            count=count+1

    init=numsave%numsectors
    n=init

    num=np.int((xarcq.size-init)/numsectors)
    n=xarcq.size-init
    for i in range(num,0,-1):

        lima=n-numsectors
        limb=n

        xradq=np.append(xradq,np.mean(xarcq[lima:limb]))
        ysbq=np.append(ysbq,np.mean(ymgeq[lima:limb]))
        ysberrq=np.append(ysberrq,np.std(ymgeq[lima:limb]))
        ysbc=np.append(ysbc,np.mean(ymgec[lima:limb]))
        ysbcerr=np.append(ysbcerr,np.std(ymgec[lima:limb]))

        n=n-numsectors

    return xradq, ysbq, ysberrq, ysbc,ysbcerr


