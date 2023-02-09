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

from ellipsect.sectors.comp import FindSB
from ellipsect.inout.plots import PlotSB

from ellipsect.sectors.num import Re90

from ellipsect.sectors.comp import SubComp


from ellipsect.inout.prt import PrintEllFilesGax
from ellipsect.inout.prt import PrintFilesGax
from ellipsect.inout.prt import PrintFilesComps

from ellipsect.inout.plots import PlotMul

from ellipsect.lib.clas import DataMge
from ellipsect.lib.clas import DataMgeComp


def EllipSectors(ellconf, galhead, galcomps, sectgalax, sectmodel, sectcomps, n_sectors = 19, minlevel = 0):


   
    #galaxy
    xradq, ysbq, ysberrq = sect2xy(sectgalax, ellconf, galhead, n_sectors)

    #model
    xradm, ysbm, ysberrm = sect2xy(sectmodel, ellconf, galhead, n_sectors)

    # Plotting
    limx,limy, axsec = PlotSB(xradq, ysbq, ysberrq, xradm, ysbm, ysberrm, ellconf, galhead.scale)

    ### surface brightness output file

    if ellconf.flagsbout == True: 

        #print to file    
        PrintEllFilesGax(ellconf, galhead, xradq, ysbq, ysberrq, xradm, ysbm, ysberrm)

    #### Creating Subcomponents images with Galfit

    if ellconf.flagcomp:

        xradq, ysbq, n = SubComp(ellconf, galhead, galcomps, sectcomps, axsec, n_sectors = n_sectors)


    axsec.legend(loc=1)

    return limx,limy


def sect2xy(sect, ellconf, galhead, n_sectors):

    #######################################
    #######################################

    # model
    stidx = np.argsort(sect.radius)

    mgerad = sect.radius[stidx]
    mgecount = sect.counts[stidx]
    mgeangle = sect.angle[stidx]
    mgeanrad = np.deg2rad(mgeangle)

    ab = ellconf.qarg

    aellab = mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)


    aellarc = aellab*galhead.scale

    # formula according to cappellary mge manual
    mgesb = galhead.mgzpt - 2.5*np.log10(mgecount/galhead.exptime) \
            + 2.5*np.log10(galhead.scale**2) + 0.1 - ellconf.Aext

    stidxq = np.argsort(aellarc)

    xarc = aellarc[stidxq]
    ymge = mgesb[stidxq]


    ######  Function to order SB along X-axis for model

    xrad, ysb, ysberr = FindSB(xarc, ymge, n_sectors)


    return xrad, ysb, ysberr




#sectors/ellipsectors.py
def MulEllipSectors(ellconf, galhead, galcomps, sectgalax, sectmodel, sectcomps):

  
    mgegal = DataMge()
    mgemod = DataMge()
    mgecom = DataMgeComp()

    eps = 1 - ellconf.qarg


    yctemp = ellconf.xc
    xctemp = ellconf.yc


    #angsec = 90 - ellconf.parg # moved to PlotMul

    ######################

    sg = sectgalax

    sm = sectmodel

    ###################################################

    stidx = np.argsort(sg.radius)

    #   galaxy
    mgegal.rad = sg.radius[stidx]

    mgegal.count = sg.counts[stidx]
    mgegal.angle = sg.angle[stidx]
    mgegal.anrad = np.deg2rad(mgegal.angle)


    # model

    stidx = np.argsort(sm.radius)

    mgemod.rad=sm.radius[stidx]

    mgemod.count=sm.counts[stidx]
    mgemod.angle=sm.angle[stidx]
    mgemod.anrad=np.deg2rad(mgemod.angle)


    # converting to pixels

    mgegal.rad = mgegal.rad*galhead.scale
    mgemod.rad = mgemod.rad*galhead.scale


    # formula according to cappellary mge manual
    # galaxy:
    mgegal.sb = galhead.mgzpt - 2.5*np.log10(mgegal.count/galhead.exptime) + 2.5*np.log10(galhead.scale**2) + 0.1

    # Model:
    mgemod.sb= galhead.mgzpt - 2.5*np.log10(mgemod.count/galhead.exptime) + 2.5*np.log10(galhead.scale**2) + 0.1


    if ellconf.flagcomp:

        angtemp=[]
        #mgesbsub=[]
        #mgeradsub=[]
        #mgeanglesub=[]
        #sectorsub=[]

        ###############################
        
        ab = ellconf.qarg
        ni=0
        maskgal = galcomps.Active == True

        while(ni < len(galcomps.N[maskgal])):

            subcmp = sectcomps[ni]


            subidx = np.argsort(subcmp.radius)

            temprad = subcmp.radius[subidx]

            #converting to arcsec
            temprad = temprad*galhead.scale

            mgecountsub=subcmp.counts[subidx]

            tempangle=subcmp.angle[subidx]
            mgeanradsub=np.deg2rad(tempangle)


            # formula according to cappellary mge manual
            tempmge= galhead.mgzpt - 2.5*np.log10(mgecountsub/galhead.exptime) + 2.5*np.log10(galhead.scale**2) + 0.1
            

            tempsectorsub = np.unique(tempangle)

            mgecom.sector.append(tempsectorsub)
            mgecom.sb.append(tempmge)
            mgecom.rad.append(temprad)
            mgecom.angle.append(tempangle)
            mgecom.anrad.append(mgeanradsub)

            ni+=1



    #plotting
    PlotMul(ellconf, galhead, galcomps, mgegal, mgemod, mgecom)


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

    num=np.int32((xarcq.size-init)/numsectors)
    n=xarcq.size-init
    for i in range(num,0,-1):

        lima=n-numsectors
        limb=n

        xradq=np.append(xradq,np.mean(xarcq[lima:limb]))
        ysbq=np.append(ysbq,np.mean(ymgeq[lima:limb]))
        #ysberrq=np.append(ysberrq,np.std(ymgeq[lima:limb]))
        ysberrq=np.append(ysberrq,stats.sem(ymgeq[lima:limb])) #standard error is more appropiated
        ysbc=np.append(ysbc,np.mean(ymgec[lima:limb]))
        #ysbcerr=np.append(ysbcerr,np.std(ymgec[lima:limb]))
        ysbcerr=np.append(ysbcerr,stats.sem(ymgec[lima:limb]))

        n=n-numsectors

    return xradq, ysbq, ysberrq, ysbc,ysbcerr


