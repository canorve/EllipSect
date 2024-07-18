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


from ellipsect.inout.prt import PrintEllFilesComps 
from ellipsect.inout.plots import PlotSub 

from ellipsect.sectors import ellip 

def SubComp(ellconf, galhead, galcomps, sectcomps, axsec, n_sectors=19):

    N=len(galcomps.N)

    #color value
    values = np.arange(N)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    ####################
    ab=ellconf.qarg
    n=0

    idxcol=0
    maskcomp = galcomps.Active == True
    while(n<N):

        if galcomps.Active[n] == True:
            namec = galcomps.NameComp[n] #check if name coincide, it must be

            scmp = sectcomps[idxcol] 


            ###################################################

            xradq, ysbq, ysberrq = ellip.sect2xy(scmp, ellconf, galhead, n_sectors)

            colorVal = scalarMap.to_rgba(values[idxcol])

            PlotSub(xradq, ysbq, n, axsec, namec, colorVal)


            if ellconf.flagsbout == True:
                ncomp=n+1
                ncomp=str(ncomp)

                PrintEllFilesComps(ellconf, galhead, namec, ncomp, xradq, ysbq, ysberrq)

            idxcol+=1
        n=n+1


    return  xradq, ysbq, n

def FindSB(xarcq, ymgeq, numsectors):
    # the xarcq array must be ordered
    # use mag instead of counts

    xradq=[]
    ysbq=[]
    ysberrq=[]
    xradq=np.array(xradq)
    ysbq=np.array(ysbq)
    ysberrq=np.array(ysberrq)

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
        ysberrq=np.append(ysberrq,stats.sem(ymgeq[lima:limb])) #standard error is more appropiate

        n=n-numsectors

    return xradq, ysbq, ysberrq


