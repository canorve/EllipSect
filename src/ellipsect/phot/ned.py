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



from ellipsect.lib.clas import DataNed

def NED(ellconf):
    "connect to NED database to obtain Gal Extinction and other variables"

    dataned = DataNed()


    objname=ellconf.objname
    band=ellconf.band

    # ignore warnings from lecture of XML file
    if not sys.warnoptions:
        warnings.simplefilter("ignore")


    ellconf.flagweb=True

    objname=ellconf.objname

    nedweb="https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname="



    if ellconf.flagnedfile:
        filened = ellconf.nedfile
    else:
        filened=ellconf.namened

    if ellconf.flagmod or ellconf.flagmag or ellconf.flagscale or ellconf.flagdim:

        GalExt=ellconf.InMagCor
        DistMod=ellconf.InDistMod
        DistMod2=ellconf.InDistMod
        Scalekpc=ellconf.InScale
        SbDim=ellconf.InSbDim
    else:
        if(not(os.path.isfile(filened))):
        #verifica si el archivo existe para no hacer conexion a internet

            if ellconf.flagnedfile:
                print("can't find user's ned {} file ".format(filened))
                ellconf.flagweb=False
            else:
                # command for wget
                #wget -O NED_51.xml "https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname=m+51"

                #https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname=ngc+6166&hconst=67.8&omegam=0.308&omegav=0.692


                wgetcmd = 'wget -O {} "{}{}&hconst={:.1f}\&omegam={:.3f}\&omegav={:.3f}"'.format(
                            filened,nedweb,objname,ellconf.hconst,ellconf.omegam,ellconf.omegav)

                print("Running: ",wgetcmd)

                print('Cosmology: hconst = {}, omega m= {}, omega lambda= {}'.format(ellconf.hconst,ellconf.omegam,ellconf.omegav)) 


                errwg = sp.run([wgetcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)

                if errwg.returncode != 0:
                    print("can't connect to NED webserver. Is your internet connection working? ")
                    print("Luminosity and absolute magnitude will not be computed") 
                    ellconf.flagweb=False
        else:
            if (ellconf.flagnedfile):
                print("using existing user's {} file ".format(filened))
            else:
                print("using existing {} file ".format(filened))
 


        print("reading ",filened)
        votable=parse(filened,pedantic=False) 

        try: 
            table=votable.get_table_by_index(0) 

        except: 
            print("I can't read file or object name can be found in file")
            print("check object name or delete NED file for a new web search")
            print("luminosity and absolute magnitude will not be computed")
            ellconf.flagweb=False 


        if ellconf.flagweb==True:
            # si flag == True
            tablephot=votable.get_table_by_id("NED_DerivedValuesTable") 

            dataphot=tablephot.array

            lumdist=dataphot["luminosity_distance"].data[0] # units in Mpc

            #DistMod= 5 * np.log10(lumdist/10) 

            DistMod=dataphot["luminosity_distance_moduli"].data[0] # units in magnitudes

            tablext=votable.get_table_by_id("NED_BasicDataTable") 
            dataext=tablext.array

            extband="gal_extinc_" + band 
            try: 
                GalExt=dataext[extband].data[0] 
            except: 
                print("can't found {} in {} check filter name. GalExt=0 ".format(extband,filened))           
                GalExt=0


            print("Luminosity distance: (Mpc) ",lumdist)

            print("Module Distance (mag): ",DistMod)


            print("Galactic Extinction for band {} : {}".format(band,GalExt))


            Scalekpc=dataphot["cosmology_corrected_scale_kpc/arcsec"].data[0] # to convert to kpc

            print("Scale kpc/arcsec",Scalekpc)


            SbDim=dataphot["surface_brightness_dimming_mag"].data[0] 

            print("SB dimming in mag",SbDim)

            ## modulo de distancia calculado en forma independiente del redshift: 
            tabledist=votable.get_table_by_id("Redshift_IndependentDistances") 

            datadist=tabledist.array

            try:
                DistMod2=datadist["DistanceModulus"].data[0] 
                DistMod2=float(DistMod2)
                print("Distance Modulus (z independent) ",DistMod2)

            except:
                print("Distance Modulus, indep. of z, can't be extracted. ")
                DistMod2=DistMod
                print("Distance Modulus will be used from luminosity distance  ",DistMod2)


        else:
            GalExt=0
            DistMod=0
            DistMod2=0
            Scalekpc=0
            SbDim=0

    # returns warnings to normal
    if not sys.warnoptions:
        warnings.simplefilter("default")

    dataned.GalExt = GalExt
    dataned.DistMod = DistMod
    dataned.DistMod2 = DistMod2
    dataned.Scalekpc = Scalekpc
    dataned.SbDim = SbDim

    return dataned 



