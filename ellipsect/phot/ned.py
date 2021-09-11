
from ellipsect.lib.libs import *

from ellipsect import *





def NED(params, galpar, galcomps):
    "connect to NED database to obtain Gal Extinction and other variables"
    
    objname=params.objname
    band=params.band

    # ignore warnings from lecture of XML file
    if not sys.warnoptions:
        warnings.simplefilter("ignore")


    params.flagweb=True

    objname=params.objname

    nedweb="https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname="

    if params.flagnedfile:
        filened = params.nedfile
    else:
        filened=params.namened

    if params.flagmod or params.flagmag or params.flagscale or params.flagdim:

        GalExt=params.InMagCor
        DistMod=params.InDistMod
        DistMod2=params.InDistMod
        Scalekpc=params.InScale
        SbDim=params.InSbDim
    else:
        #checar si el archivo existe para no hacer conexion a internet
        if(not(os.path.isfile(filened))):

            if params.flagnedfile:
                print("can't find user's ned {} file ".format(filened))
                params.flagweb=False
            else:
                # command for wget
                #wget -O NED_51.xml "https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname=m+51"
                wgetcmd = 'wget -O {} "{}{}"'.format(filened,nedweb,objname)

                print("Running: ",wgetcmd)

                errwg = sp.run([wgetcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)

                if errwg.returncode != 0:
                    print("can't connect to NED webserver. Is your internet connection working? ")
                    print("Luminosity and absolute magnitude will not be computed") 
                    params.flagweb=False
        else:
            if (params.flagnedfile):
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
            params.flagweb=False 


        if params.flagweb==True:
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



    return (GalExt,DistMod,DistMod2,Scalekpc,SbDim)



