
from ellipsect.lib.libs import *


from ellipsect.inout.prt import PrintEllFilesComps 
from ellipsect.inout.plots import PlotSub 

from ellipsect.sectors import ellip 

def SubComp(ellconf, galhead, galcomps, sectcomps, axsec, n_sectors=19):

    N=len(galcomps.N)

    #color value
    values = range(N)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    ####################
    ab=ellconf.qarg
    n=0

    maskcomp = galcomps.Active == True
    while(n<N):

        namec = galcomps.NameComp[maskcomp][n] #check if name coincide, it must be

        scmp = sectcomps[n] 


        ###################################################

        xradq, ysbq, ysberrq = ellip.sect2xy(scmp, ellconf, galhead, n_sectors)

        colorVal = scalarMap.to_rgba(values[n])

        PlotSub(xradq, ysbq, n, axsec, namec, colorVal)


        if ellconf.flagsbout == True:
            ncomp=n+1
            ncomp=str(ncomp)

            PrintEllFilesComps(ellconf, galhead, namec, ncomp, xradq, ysbq, ysberrq)


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


