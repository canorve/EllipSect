___

# **EllipSect**

[![DOI](https://zenodo.org/badge/282223217.svg)](https://zenodo.org/badge/latestdoi/282223217)

EllipSect.py creates surface brightness profiles and extracts 
other photometric data from the 
from the GALFIT output peng et al. (2002). 

___

## **Code**:

**[EllipSect.py](EllipSect.py)**

This code is "similar" (but not sustitute) to IRAF's ellipse routine. It 
creates a Surface brightness profile for the galaxy and the model.

In addition, *EllipSect* can compute variables such as Absolute Magnitude, 
luminosity, Flux, total apparent magnitude, Bulge to Total Ratio, Tidal, Chinu
in the sectors ellipse, Bumpiness, Signal to Noise Ratio, Akaike Information criterion, 
Bayesian information criterion, mean surface brightness at effective radius, percentage 
of total light per component, radius at 90% of light (for Sersic component only), 
effective radius in kpc, etc.  

___

## **Installation**

The code is written for python 3.

The python libraries used are:

- numpy
- astropy
- scipy
- matplotlib
- mgefit

Download the release and installed via

```
cd ellipsect
pip install . 
```

or 

```
cd ellipsect
python setup.py install
```

Run the automated tests:

```
python -m pytest 
```



**EllipSect** needs the GALFIT output files (GALFIT.XX) to work.
Although **GALFIT** is not stricly required, it will need it 
to create the GALFIT components and sigma image. Make it
sure that you can call GALFIT from the command line.


EllipSect uses the mgefit library which 
is described in Cappellari, MNRAS, 333, 400 (2002).

**Install the mgefit library via pip:**  

```
pip install mgefit
```

___

## **Basic run:**

Once installed, run the *ellsec.py* in the same directory 
that you run GALFIT. It only requires the latest GALFIT's 
output file. *The easiest way to run the program is:*

```
./ellsec.py galfit.01
```

It will display images like the ones below:

   ![A85 ](img/A85.def.png)

___

For linux, you can make a simbolic link to *ellsec.py*

```
cd /usr/local/bin
sudo ln -s /path/to/your/installed/ellipsect/ellsec.py .
```

for more options:

```
./ellsec.py -help 
```


### **HOW TO USE**

To see other forms of how to run and other options available see:

   [How to use](docs/howto.md)

___

## **Questions?**

Do you have any questions or suggestions?
Please send an email to canorve [at] gmail [dot] com 
or open an [issue](https://github.com/canorve/EllipSect/issues)

I'm open to new ideas that can benefit the 
software *EllipSect* and the *GALFIT* community

___

## **License**

The code is under the license of **GNU**

___

## **References**


Cappellari, MNRAS, 333, 400 (2002).


___

Check my others GALFIT tools [here](https://github.com/canorve/GALFITools)

___
