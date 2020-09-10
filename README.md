
# EllipSect

EllipSect.py creates surface brightness profiles
for galaxy and galfit model from the galfit output: galfit.XX. 
See peng et al. (2002). It also creates multiple surface brightness 
plots for separate angles.  

## Installation

Download  or clone this code. These codes are
written for python 3.

The python libraries used are:
- numpy
- sys
- os
- subprocess
- astropy
- scipy
- matplotlib
- mgefit
- mimetypes
- warnings
- platform

Although GALFIT is not *strictly* required, EllipSect needs the GALFIT output files (GALFIT.XX) to work.


## Code:

**[EllipSect.py](EllipSect.py)**

This code is similar to IRAF's ellipse routine. It 
creates a Surface brightness profile for the galaxy and model.

In addition, *EllipSect* can compute variables such as Absolute Magnitude, 
luminosity, Flux, total apparent magnitude, Bulge to Total Ratio, Tidal, Chinu
in the sectors ellipse, Bumpiness, Signal to Noise Ratio, Akaike Information criterion, 
Bayesian information criterion, mean surface brightness at effective radius, percentage 
of total light per component, radius at 90% of light (Sersic component only), 
effective radius in kpc, etc.  

## Additional Libraries: mgefit

**Install the mgefit library via pip:**  

```
pip install mgefit
```

*The easiest way to run the program is:*

```
 ./EllipSect.py galfit.01 
```

## OPTIONS

The options to run the code in the terminal (or ipython) are:

```
 ./EllipSect.py [GALFITOutputFile] [-logx] [-q AxisRatio] [-pa PositionAngle] [-comp] [-pix] [-ranx/y Value] [-grid] [-dpi Value] [-noplot] [-phot] [-sbout] [-noplot] [-minlevel Value] [-sectors Value] [-object Name] [-filter Name] [-snr] [-help] [-checkimg] [-noned] [-distmod Value] [-magcor Value] [-scalekpc Value][-sbdim Value] [-model ModelImage] [-keep] [-ned XmlFile] 
 ```

### Input File

**GALFITOutputFile**: GALFIT output file  (e.g. galfit.01)


**help**: Help menu

**logx**: plots X-axis as logarithm

**q**: axis ratio value. If ignored, it takes the one from the last component in GALFITOutputFile.

**pa**: position angle value (same as GALFIT). If ignored, it takes the one 
from the last component in GALFITOutputFile.

**comp**: plots include the individual model components

**pix**: adds pixels units to the top of x-axis.

**ranx**: constant that increase/decrease the range of the x axis *or* 
it can be used as xmin-xmax to change range.

**rany**: constant that increase/decrease the range of the y axis *or* 
it can be used as ymin-ymax to change range.

**noplot**: the code creates images but do not display them.

**grid**: display a grid in the plot

**dpi**: dots per inch value to increase/decrease resolution.

**sbout**: Creates output file containing the surface brightness profiles.


### Photometry output options

**phot**: Compute photometry. Check the variables created in output file.

The below options are used only if 'phot' is enabled:

**snr**: Creates a signal to noise image. This is created dividing the galaxy image
with the sigma image created by GALFIT

**object**: used for 'phot' to search in NED. For instance, if you are looking
for photometry data for galaxy m51, then used as "-object m51" the same name 
that you will used to search in NED.

**filter**: used for 'phot' to indicate band for NED. If you need galactic correction
for B  filter then used as "-filter B". Band "R" is the default option.  

Any of the following options disabled the connection to NED

**noned**: it avoids to connect with NED. No luminosity nor absolute magnitude is computed.

**distmod**: manual input for Distance Modulus.  

**magcor**: manual input for Galactic Extinction. 

**scalekpc**: manual input for equivalence of ''/kiloparsec. 

**sbdim**: manual input for surface brightness dimming.

### Advanced 

**minlevel**: Parameter given directly to sectors_photometry.
              Ellipse radius stops when it founds this value. Check sectors_photometry manual 

**sectors**: parameter given directly to sectors_photometry. Divide ellipse in 'sectors'
                      Check sectors_photometry manual
                     
**checkimg**: save the images used in sectors_photometry for individual components

**model**: User can introduce his/her own model image. SNR quantities will be 
            inaccurate  for  this option.

**sky**: User can introduce his/her own sky value. 

**keep**: Use existing '-comp.fits' file to compute subcomponents. 

**ned**: User can introduce his/her own  ned xml file.



## Notes

* EllipSect works for GALFIT version >*3.0.7* Check this because you will
    have issues for *-comp* option.

* EllipSect uses the mask image (option "*F*" GALFIT) only if this
   is a **FITS** image. In case your mask is an *ASCII* file, you can convert it to **FITS** using the [xy2fits.py](https://github.com/canorve/GALFITools/blob/master/docs/xy2fits.md) tool.


* EllipSect uses axis ratio (*q*) and position angle (*pa*) to create an "ellipse" *grid* using the function *sectors_photometry* from the *mgefit* library. Unlike IRAF's Ellipse, *q* and *pa* are fixed through radius. See the images below: 

    ![A85 gal](img/A85.png)

    For this reason, errors are expected to 
    be greater than those coming from IRAF's ellipse since EllipSect 
    averages errors for different isophotes. While, on the other hand, IRAF's ellipse
    can change axis ratio and angular position for each isophote. 
    
    This is how mgefit *sectors_photometry* returns the counts data and, unless I write my own code, I can't change that. 

* Be sure to run this code in the same path that you run GALFIT. 

* Since EllipSect reads the "*B*" option of galfit.XX file, this
    must be the **last** GALFIT fit. 

* The angles shown in the multi-plot are measured from the galaxy's major axis.
    They are **not** measured from the Y-axis. 

* In order for the program to detect the components, they must share the same 
    center (x,y). This allows that galfit can use other components such as
     sersic to be used as masks for nearby GALFIT. EllipSect already does not 
     take them into account for the plots.

* For the comp option, It could be some small differences between the angle shown 
    in the top right corner and the one from each component. This is because
    *sectors_photometry* is applied different for individual components and the 
    galaxy itself. They are at different angles. To see the real angle which the 
    component is measured check the output file at that angle with the *-sbout* option 

## Examples

* Displays the help menu: 
    ./EllipSect.py -help


* To introduce an axis ratio of 0.35 and position angular of 60 
    (measured from Y-axis): 

     ./EllipSect.py galfit.02 -q 0.35 -pa 60 

### Plot Examples

See the examples below for an elliptical galaxy that was fitted 
with 7 gaussians (images for this galaxy are displayed above). 

* Simple plot example: 
    ./EllipSect.py galfit.46 

    ![A85 ](img/A85.def.png)
 
* Displays the X-axis as log:
    ./EllipSect.py galfit.46 -logx

    ![A85 ](img/A85.log.png)

* Include the individual model components to the plot:  
    ./EllipSect.py galfit.46 -comp
    (displays the 7 gaussians)

    ![A85 ](img/A85.comp.png)

* Insert pixels units in the top X-axis: 
    ./EllipSect.py galfit.46 -pix

    ![A85 ](img/A85.pix.png)

* Range in X-axis is decreased 50%: 
    ./EllipSect.py galfit.46 -ranx 0.5 

    ![A85 ](img/A85.ranx1.png)

* Range in Y-axis is doubled: 
    ./EllipSect.py galfit.46 -rany 2 

    ![A85 ](img/A85.rany1.png)

* X-axis range vary from 1 to 50: 
    ./EllipSect.py galfit.46 -ranx 1-50 

    ![A85 ](img/A85.ranx2.png)

* Use grid on plot and increase resolution to 300 dots per inch: 
    ./EllipSect.py galfit.46 -grid -dpi 300 

    ![A85 ](img/A85.grid.png)

* Same as above but popup window does not appear. Plots are 
    directly saved in directory: 
 
    ./EllipSect.py galfit.46 -grid -dpi 300 -noplot 

* If the user desires to create their own plots, 'sbout' option
  will save the surface brightness vs. radius of the galaxy and model
  in a file:
   
    ./EllipSect.py galfit.46 -sbout

  EllipSect can also save the surface brightness data for 
  individual components in separated files:

    ./EllipSect.py galfit.46 -comp -sbout 

### Phot Examples

*   *EllipSect* can calculate additional info besides the ones 
    that are already included in the galfit.XX or fit.log files.  
    Those variables are intended to help the user to have a quick reference 
    of the model and decide to modify the model or increase the number of components.
    
    Such output variables have to be taken with caution and they always have to be verified for the user before to include them in the final paper.

    The output photometry variables include: Absolute Magnitude, luminosity, Flux, total apparent magnitude, Bulge to Total Ratio, Tidal, $\Chi_\nu$ within the sectors ellipse, Bumpiness, Signal to Noise Ratio, Akaike Information Criterion, Bayesian Information Criterion, mean surface brightness, percentage of total light per individual component, radius at 90% of light ( for Sersic components only). 
    
    Those variable are stored in a single file when the following command is executed:

    ./EllipSect.py galfit.46 -phot


*   *phot* option looks in NED (NASA/IPAC Extragalactic Database) for 
    Galactic Extinction, distance modulus, surface brightness dimming, etc. 
    to compute Absolute Magnitude, luminosity and other variables. To do this,
    EllipSect looks for name of the galaxy (as it is searched in NED) 
    and wavelength band in the header. If that info is not in 
    the header, the user can introduce the band and object name as
    it is shown in the next example for galaxy messier 51 in the band B:  


    ./EllipSect.py galfit.14 -phot -object m51 -filter B

*   If the user wants to see a Signal to Noise image of the data, use 
    the next command:

    ./EllipSect.py galfit.14 -phot -snr

*   If for some reason the user does not want to connect to NED use 
    the following option:

    ./EllipSect.py galfit.14 -phot -noned

    take into account that Luminosity and Absolute magnitud will not be computed

*   EllipSect allows to introduce manually the NED info. For example, 
    the next command introduce a distance modulus of 10, galactic extinction 
    of 0.3, "/kpc of 1.3 and surface brightness dimming of 0.3.  

    ./EllipSect.py galfit.10 -phot -distmod 10 -magcor 0.3 -scalekpc 1.3 -sbdim .3

    This option avoids to connect with NED. 

    EllipSect does not correct by K-correction. 


### Advanced Examples

*   model option allows the user to introduce his own model image for analysis.
    EllipSect will use this image instead of the one created by GALFIT output.
    If -phot option is enabled, SNR quantities will be inaccurate.

    ./EllipSect.py galfit.14 -model model.fits

*   sky option allows the user to introduce his/her own sky value to substract 
    from galaxy and model images. User must know what he/she is doing because
    otherwise EllipSect will produce wrong outputs.  

    ./EllipSect.py galfit.01 -sky 300


*   The following options requires that the user has already experienced with 
    the *sectors_photometry* function of the mge library. 

    minlevel is a parameter that is given directly to *sectors_photometry* 
    It indicates when the functions stops. For example, the following command
    tells to *sectors_photometry* that stops when the sky is 0.

    ./EllipSect.py galfit.14 -minlevel 0


    Note: Galfit sky parameter is already removed from image before the call 
    to *sectors_photometry* 


*   sectors option is another parameter that is given directly to *sectors_photometry*.
    It tells the function in how many sectors it should divide. *sectors_photometry* 
    use four-fold symmetry.  

    ./EllipSect.py galfit.14 -sectors 19

*   checkimg will create images used by *sectors_photometry* to check how 
    this function was used on the individual model components. The images 
    names will start with 'C' followed by the component number. 
    
    Use it with the 'comp' option:  

    ./EllipSect.py galfit.14 -comp -checkimg

### Run the example

Check the [example](example) to test it by yourself.








**EllipSect.py** uses the mgefit library which is
described in Cappellari, MNRAS, 333, 400 (2002).










## Questions?

Do you have any questions or suggestions?
Please send an email to canorve [at] gmail [dot] com 
or open an [issue](https://github.com/canorve/EllipSect/issues)

I'm open to new ideas that can benefit the GALFIT community

## License

These codes are under the license of **GNU**

