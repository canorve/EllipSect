from ellipsect.lib.libs import *

from ellipsect import *



def Help():

    print ("EllipSect: an analysis tool for GALFIT output  ")
    print ("Usage:\n %s [GALFITOutputFile] [-logx] [-q AxisRatio] [-pa PositionAngle] [-comp] [-pix] [-ranx/y Value] [-grid] [-dpi Value] [-model File] [-phot] [-sky Value] " % (sys.argv[0]))
    print ("More options: [-sbout] [-noplot] [-minlevel Value] [-sectors Value] [-object Name] [-filter Name] [-snr] [-help] [-checkimg] [-noned] [-distmod Value] [-magcor Value] [-scalekpc Value] [-sbdim Value] [-keep] [-ned XmlFile] [-gradsky ] [-randsky ] [-skyRad Value] [-skyRadmax Value][-skynum Value] [-skybox Value] [-skywidth Value] [-fwhm Value]") 

    print ("GALFITOutputFile: GALFIT output file ")
    print ("logx: activates X-axis as logarithm ")
    print ("q: introduce axis ratio ")
    print ("pa: introduce position angle (same as GALFIT) ")
    print ("comp: plots individual components ")
    print ("pix: plot the top x-axis in pixels ")
    print ("ranx: constant that increases the range of the x axis (for arcsec only) ")
    print ("another way to use ranx is to provide a range: xmin-xmax")
    print ("rany: constant that increases the range of the y axis (surface brightness)")
    print ("another way to use ranx is to provide a range: ymin-ymax")
    print ("grid: display a grid in the plot ")
    print ("dpi: dots per inch used for images files ")
    print ("noplot: avoid displaying windows and directly creates images")
    print ("sbout: creates output file containing the surface brightness profiles")
    print ("       All surface brightness files will be saved in 'sbfiles' directory")
    print ("keep: use existing file to compute subcomponents ")
        

    print ("                OUTPUT               ")
    print ("phot: Compute photometry. Check the created output file")
    print ("the below options are used only if 'phot' is enabled ")    
    print ("      snr: Creates Signal to Noise image ")
    print ("      object: used for 'phot' to search in NED  ")
    print ("      filter: used for 'phot' to indicate band for NED ")
    print ("      ned: user can introduce his/her own ned xml file")

    print ("      noned: avoid to connect to NED")
    print ("any of the following options disabled the connection of NED ")    
    print ("      distmod: Introduce Distance Modulus ")
    print ("      magcor: Introduce Galactic Extinction ")
    print ("      scalekpc: Introduce equivalence of ''/kiloparsec ")
    print ("      sbdim: Introduce surface brightness dimming")
    print ("                ADVANCED               ")
    print ("model: User can introduce his/her own image model.")
    print ("sky: User can introduce his/her own sky value.")
    print ("minlevel: parameter given directly to sectors_photometry.")
    print ("                      It stops when it founds this value")

    print ("sectors: parameter given directly to sectors_photometry. Divide elipse in 'sectors' ")
    print ("                      Check sectors_photometry manual")
    print ("checkimg: save the images used for sectors_photometry in individual components")
    print ("fwhm: It is used to compute Area_psf for BICres. Default = 2 pixels")

    print ("gradsky: computes sky using the gradient method ")
    print ("randsky: computes sky averaging random boxes ") 

    print ("skyRad: for randsky, it creates a mask for the main target using this radio. For gradsky it is where the program starts to compute the gradient.")
    print ("skyRadmax: for randsky only, maximum radius from main target where randbox can be selected ")
    print ("skynum: Number of boxes used in randsky. Default = 20")
    print ("skywidth: width of the ring for gradsky. Default = 20")
    print ("skybox: pixel size of the box for randsky. Default = 20") 

 
    print ("\n ")
    print ("help: This menu ")
    print ("\n ")

    print ("Example:\n %s galfit.01 -logx" % (sys.argv[0]))
    print ("or Example:\n %s galfit.02 -q 0.35 -pa 60 -comp -ranx 2 -phot " % (sys.argv[0]))
    print ("or Example:\n %s galfit.02 -q 0.35 -pa 60 -comp -ranx 1-20 \n" % (sys.argv[0]))
    print ("see https://github.com/canorve/EllipSect/tree/master/docs/howto.md   for more examples")



    sys.exit()


    return True



