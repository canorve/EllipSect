
from ellipsect.lib.libs import *

from ellipsect import *



class SkyCal:

    def RandBox(self,ImageFile,MaskFile,xx,yy,thetadeg,q,Rinit,box,num,Rmax):

        self.xx = xx 
        self.yy = yy

        #angsec=90-galpar.ang
        #self.thetadeg= thetadeg
        #self.thetadeg = 90 - thetadeg
        self.thetadeg=90 + thetadeg

        self.q = q
        self.Rinit = Rinit
       
        self.box = box 
        self.num = num

        ###

        hdumask = fits.open(MaskFile)
        self.maskimg = hdumask[0].data
        hdumask.close()


        hdu=fits.open(ImageFile)
        self.img = hdu[0].data
        hdu.close()

        ####
        
        (self.nrow,self.ncol)=self.img.shape

        xmin,xmax,ymin,ymax,Rend = self.GetXYRBorder()

        if Rmax != 0:
            if Rmax <= Rend:
                Rend = Rmax
            else:
                print("input skyRadmax is greater than image size")
            print("using Rmax = ",Rend)

        #print("border x ",xmin,xmax)
        #print("border y ",ymin,ymax)
        #print("Rborder ",Rkron)

        mean, std, median  = self.GetRandBoxSky( Rinit, Rend )

        #print("Total sky:  mean, std, median = ",mean,std,median)

        meansky = np.mean(mean)
        medsky = np.median(median)
        rmstd = np.sqrt(np.mean(std**2))



        return meansky, rmstd, medsky


    def GetXYRBorder(self):
        "this subroutine get the coordinates of the border"

        q =  self.q

        theta = self.thetadeg * (np.pi / 180)  # rads!!

        thetax=np.sqrt((np.cos(theta))**2 + (q**2)*(np.sin(theta))**2 )
        thetay=np.sqrt((q**2)*(np.cos(theta))**2 + (np.sin(theta))**2 )


        if (self.thetadeg >-45 and self.thetadeg <= 45):

            xmax=self.ncol
            xmin =1
            R1 = (xmax - self.xx)/thetax
            R2 = (self.xx - xmin)/thetax

            if (np.abs(R1)>=np.abs(R2)):
                R = R1
            else:
                R = R2

        elif (self.thetadeg >45 and self.thetadeg <= 135):
            ymax=self.nrow
            ymin =1
            R1 = (ymax - self.yy)/thetay
            R2 = (self.yy - ymin)/thetay

            if (np.abs(R1)>=np.abs(R2)):
                R = R1
            else:
                R = R2

        elif (self.thetadeg >135 and self.thetadeg <= 225):
            xmax=1
            xmin =self.ncol

            R1 = (xmax - self.xx)/thetax
            R2 = (self.xx - xmin)/thetax

            if (np.abs(R1)>=np.abs(R2)):
                R = R1
            else:
                R = R2

     
        elif (self.thetadeg >225 and self.thetadeg <= 315):
            ymax=1
            ymin =self.nrow
            R1 = (ymax - self.yy)/thetay
            R2 = (self.yy - ymin)/thetay

            if (np.abs(R1)>=np.abs(R2)):
               R = R1
            else:
               R = R2

        R = np.abs(R) #avoids negative numbers
        bim = q * R
        # getting size

        xmin = self.xx - np.sqrt((R**2) * (np.cos(theta))**2 +
                         (bim**2) * (np.sin(theta))**2)

        xmax = self.xx + np.sqrt((R**2) * (np.cos(theta))**2 +
                         (bim**2) * (np.sin(theta))**2)

        ymin = self.yy - np.sqrt((R**2) * (np.sin(theta))**2 +
                         (bim**2) * (np.cos(theta))**2)

        ymax = self.yy + np.sqrt((R**2) * (np.sin(theta))**2 +
                         (bim**2) * (np.cos(theta))**2)

        mask = xmin < 1
        if mask.any():
            if isinstance(xmin,np.ndarray):
                xmin[mask] = 1
            else:
                xmin = 1

        mask = xmax > self.ncol
        if mask.any():
            if isinstance(xmax,np.ndarray):
                xmax[mask] = self.ncol
            else:
                xmax = self.ncol

        mask = ymin < 1
        if mask.any():
            if isinstance(ymin,np.ndarray):
                ymin[mask] = 1
            else:
                ymin = 1

        mask = ymax > self.nrow
        if mask.any():
            if isinstance(ymax,np.ndarray):
                ymax[mask] = self.nrow
            else:
                ymax =self.nrow

        xmin=np.int(np.round(xmin))
        ymin=np.int(np.round(ymin))
        xmax=np.int(np.round(xmax))
        ymax=np.int(np.round(ymax))


        return (xmin,xmax,ymin,ymax,np.int(R))


    def GetRandBoxSky(self, Rinit, Rmax):

        #print("Rinit, Rmax ",Rinit, Rmax)
        (nrow,ncol)=self.img.shape
        # it obtains corners of Rinit
        (xmino, xmaxo, ymino, ymaxo) = self.GetSize(self.xx, self.yy, Rinit, self.thetadeg, self.q, self.ncol, self.nrow) 
        # it obtains corners of Rmax
        (xminf, xmaxf, yminf, ymaxf) = self.GetSize(self.xx, self.yy, Rmax, self.thetadeg, self.q, self.ncol, self.nrow) 
        
        Value=1 #  value of counts  of the  main target for  the mask image  
        self.maskimg = self.MakeKron(self.maskimg, Value, self.xx, self.yy, Rinit, self.thetadeg, self.q, xminf, xmaxf, yminf, ymaxf) 

        ########

        sky = np.array([])
        skystd = np.array([])
        skymed = np.array([])


        cont=self.num # este es el numero de cajas que va a utilizar 
        for idx,item in enumerate(range(cont)):

            flatimg,xinit,yinit=self.GetRandomPatch(self.img,self.maskimg,self.box,xmino,xmaxo,ymino,ymaxo,xminf,xmaxf,yminf,ymaxf)
            xfin = xinit + self.box - 1
            yfin = yinit + self.box - 1 


            flatimg.sort()

            boxcont=0

            while( not(flatimg.any()) and (boxcont < 10)):

                if (boxcont == 0):
                    print("Picking another box ")

                flatimg,xinit,yinit=self.GetRandomPatch(self.img,self.maskimg,self.box,xmino,xmaxo,ymino,ymaxo,xminf,xmaxf,yminf,ymaxf)
                xfin = xinit + self.box - 1
                yfin = yinit + self.box - 1


                flatimg.sort()

                boxcont+=1 # avoid eternal loop

            if (boxcont == 10):
                print("max. iteration reached. I couldn't find a box")
     
            tot=len(flatimg)

            top=round(.8*tot)
            bot=round(.2*tot)
            
            imgpatch=flatimg[bot:top]

            linebox = "Box:{}   xinit/fin: {}-{}; yinit/fin: {}-{}  ".format(item+1,xinit,xfin,yinit,yfin)
            #print(linebox)

            mean=np.mean(imgpatch)
            std=np.std(imgpatch)

            median=np.median(imgpatch)

            linemean = "sky  mean = {:.2f}; std = {:.2f}; median = {:.2f}".format(mean,std,median)
            print(linemean)

            sky=np.append(sky,mean)
            skystd=np.append(skystd,std)
            skymed=np.append(skymed,median)

        return sky,skystd,skymed


    def MakeKron(self,imagemat, idn, x, y, R, theta, q, xmin, xmax, ymin, ymax):
        "This subroutine create a Kron ellipse within a box defined by: xmin, xmax, ymin, ymax"

        xmin = np.int(xmin)
        xmax = np.int(xmax)
        ymin = np.int(ymin)
        ymax = np.int(ymax)

        #q = (1 - ell)
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

    def GetSize(self,x, y, R, theta, q, ncol, nrow):
        "this subroutine get the maximun"
        "and minimim pixels for Kron and sky ellipse"
        # k Check
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

        xmin=np.int(np.round(xmin))
        ymin=np.int(np.round(ymin))
        xmax=np.int(np.round(xmax))
        ymax=np.int(np.round(ymax))


        return (xmin, xmax, ymin, ymax)

    def GetRandomPatch(self,imagemat,mimg,box,xmino,xmaxo,ymino,ymaxo,xminf,xmaxf,yminf,ymaxf):

        # get a random box patch of the imagemat 
        xinit,yinit=self.GetRandomCoord(xmino-box,xmaxo+box,ymino-box,ymaxo+box,xminf,xmaxf,yminf,ymaxf)

        xfin = xinit + box  - 1 
        yfin = yinit + box  - 1

        imagebox=imagemat[yinit - 1:yfin, xinit - 1:xfin]
        maskbox=mimg[yinit - 1:yfin, xinit - 1:xfin]

        invboxpatch=np.logical_not(maskbox)

        return imagebox[invboxpatch],xinit,yinit



    def GetRandomCoord(self,xmino,xmaxo,ymino,ymaxo,xminf,xmaxf,yminf,ymaxf):
        #obtains coordinates between the inner and outer box
        coordinates = [(x,y) for x in np.arange(0,xmaxf) for y in np.arange(0,ymaxf) if not((x >= xmino and x <= xmaxo) and ( y >= ymino and y <= ymaxo))]

        # choose xinit, and yinit from coordinates

        ridx = np.random.randint(0,len(coordinates)-1)

        xinit = coordinates[ridx][0] 
        yinit = coordinates[ridx][1] 

        return xinit,yinit 

        ######


    def GetEllipSky(self, ImageFile, MaskFile, xx, yy, thetadeg, q, Rinit, width,namering):

        self.xx = xx 
        self.yy = yy

        #angsec=90-galpar.ang
        #self.thetadeg= thetadeg
        #Theta=galpar.ang + 90
        self.thetadeg=90 + thetadeg
        self.q = q
        self.e = (1 - self.q)
        self.Rinit=Rinit
       
        self.width = width 

        # check how many galapagos has to compute sky gradient num = 15?
        self.NumRings = 5  # number of rings per loop # read this from function?
        ###


        hdumask = fits.open(MaskFile)
        self.maskimg = hdumask[0].data
        hdumask.close()


        hdu=fits.open(ImageFile)
        self.img = hdu[0].data
        hdu.close()

        ####
 
       
        (self.nrow,self.ncol)=self.img.shape


        xmin,xmax,ymin,ymax,Rkron = self.GetXYRBorder()

        self.R=Rkron

        #Rinit=2 # init in 2 pixels 


        #bim = q * self.R   # check this bim and below, is R or Rinit
        #bim = q * Rkron   # check this bim and below, is R or Rinit


        #R2 = self.R + self.width # not needed?
        #bim2 = q * R2


#        (xmino, xmaxo, ymino, ymaxo) = self.GetSize(self.xx, self.yy, Rinit, self.thetadeg, self.e, self.ncol, self.nrow) 

        (xmin, xmax, ymin, ymax) = self.GetSize(self.xx, self.yy, Rkron, self.thetadeg, self.q, self.ncol, self.nrow) # obtain corners of R

        theta = self.thetadeg * np.pi / 180  # Rads!!!

        ypos, xpos = np.mgrid[ymin - 1:ymax, xmin - 1:xmax]


        patch = self.maskimg[ymin - 1:ymax, xmin - 1:xmax] # logical patch mask image

        self.invpatch=np.logical_not(patch)

        dx = xpos - self.xx
        dy = ypos - self.yy

        self.dist = np.sqrt(dx**2 + dy**2)

        landa = np.arctan2(dy, dx)


        mask = landa < 0
        if mask.any():
            landa[mask] = landa[mask] + 2 * np.pi

        landa = landa - theta


        Rings=np.arange(self.Rinit,self.R,self.width) # anillos de tamaño width
        #bim=np.arange(self.Rinit*q, self.Rinit*q+self.width*len(Rings), self.width) # asi esta para que de el mismo tamanio que Rings
       
        #bRings=Rings*q


        sky = np.array([])
        skymed = np.array([])
        skystd = np.array([])
        radius = np.array([])

        flagfirst=True

        count = 0
        idx=0
        idx2 = idx + 1


        for ind, item in enumerate(Rings):


            maskring,flagfirst=self.GetRingMask(Rings,landa,theta,flagfirst, idx,idx2)


            flatimg=self.img[ypos[maskring], xpos[maskring]].flatten()  
            flatimg.sort()

            tot=len(flatimg)

            top=round(.8*tot)
            bot=round(.2*tot)


            imgpatch=flatimg[bot:top]

            mean=np.mean(imgpatch)
            std=np.std(imgpatch)

            median=np.median(imgpatch)

            sky=np.append(sky,mean)
            skymed=np.append(skymed,median)
            skystd=np.append(skystd,std)
            radius=np.append(radius,Rings[idx] + self.width/2)

            print("rad = {:.2f}; sky mean = {:.2f}; sky std = {:.2f}; median: {:.2f} ".format(Rings[idx] + self.width/2,mean,std, median))

            ###

            # calcular gradiente
            if (count >= self.NumRings):
                gradmask = np.gradient(sky[1:-1]) >= 0 # [1:-1] avoiding the first and last element for gradient 
               
                count = 0
                tempidx=np.where(np.gradient(sky[1:-1]) >= 0)
                if (sky[1:-1][gradmask].any()): 
                    
                    savidx=tempidx[0][0]
                    savidx2=savidx+1                
                    
                    flagfirst=True
                    maskring,flagfirst=self.GetRingMask(Rings[1:-1],landa,theta,flagfirst, savidx, savidx2)
                    
                    print("Ring radius = {:.2f} marked in {} ".format(radius[1:-1][savidx],namering))
                    print("In this file, the count value within ring represent the value of the long axis") 
                    self.img[ypos[maskring], xpos[maskring]] = radius[1:-1][savidx] 
                    break

                #clean the arrays:
                #sky = np.array([])
                #skymed = np.array([])
                #skystd = np.array([])
                #radius = np.array([])

            count += 1
            idx +=1
            idx2 +=1

            if idx == (len(Rings)-1): 
                print("The edge of image has been reached. Sky can not be computed")
                return 0,0,0,0


        hdu[0].data=self.img

        hdu.writeto(namering,overwrite=True) 

        finmean,finmedian,finstd,finRad = sky[1:-1][gradmask],skymed[1:-1][gradmask],skystd[1:-1][gradmask],radius[1:-1][gradmask]

        #rms = np.sqrt(np.mean(mean**2))
        #rmstd = np.sqrt(np.mean(std**2))
        #rmsmed = np.sqrt(np.mean(median**2))


        return finmean[0],finstd[0],finmedian[0],finRad[0]


    def GetRingMask(self,Rings,landa,theta,flagfirst, idx,idx2):


        bim=Rings[idx]*self.q
        tempangle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / Rings[idx])

        tempxell = self.xx + Rings[idx] * np.cos(tempangle) * np.cos(theta) - bim * \
          np.sin(tempangle) * np.sin(theta)

        tempyell = self.yy + Rings[idx] * np.cos(tempangle) * np.sin(theta) + bim * \
          np.sin(tempangle) * np.cos(theta)

        
        bim2=Rings[idx2]*self.q
        tempangle2 = np.arctan2(np.sin(landa) / bim2, np.cos(landa) / Rings[idx2])

        tempxell2 = self.xx + Rings[idx2] * np.cos(tempangle2) * np.cos(theta) - bim2 * \
          np.sin(tempangle2) * np.sin(theta)

        tempyell2 = self.yy + Rings[idx2] * np.cos(tempangle2) * np.sin(theta) + bim2 * \
          np.sin(tempangle2) * np.cos(theta)
 

        if (flagfirst):
            self.dell = np.sqrt((tempxell - self.xx)**2 + (tempyell- self.yy)**2)
            flagfirst=False
        else:
            self.dell = self.dell2

        self.dell2 = np.sqrt((tempxell2 - self.xx)**2 + (tempyell2 - self.yy)**2)
        maskring = (self.dist < self.dell2) & (self.dist > self.dell)  
        maskring=maskring*self.invpatch

        return maskring,flagfirst 


################################################
################################################
################################################


