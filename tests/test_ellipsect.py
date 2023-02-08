

import pytest
import os


import subprocess as sp


from ellipsect import ellipsectors 

from ellipsect.inout.read import ArgParsing
from ellipsect.sectors.sect import SectorsGalfit


# simple run test
def test_exit():
    with pytest.raises(SystemExit) as e:
        ellipsectors.run()
    assert e.type == SystemExit 
    assert e.value.code == 2 


#check if galfit is installed
def test_galfit():

    runcmd="galfit -help"
    err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  

    assert err.returncode == 0, "is GALFIT installed?"




# checking the creation of files
def test_files():

    arg=['tests/galfit.01', '-np']

    path="tests/"

    filepng = "imgblock.png"
    filemulpng = "imgblock-mul.png"

    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)



    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)



    assert os.path.isfile(filepng)
    assert os.path.isfile(filemulpng)


    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)




# checking the creation of the components file
def test_comp():

    arg=['tests/galfit.01','--comp', '--noplot']

    path="tests/"

    filepng = "imgblock.png"
    filemulpng = "imgblock-mul.png"

    filepng = path+filepng
    filemulpng = path+filemulpng



    filecomp = "imgblock-comp.fits"

    filecomp= path+filecomp

    if os.path.isfile(filecomp):
        os.remove(filecomp)


    # read user's input 
    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)




    assert os.path.isfile(filecomp),"is GALFIT installed?"


    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)


    if os.path.isfile(filecomp):
        os.remove(filecomp)






# checking the creation of the sbout files
def test_phot():


    arg=['tests/galfit.01','--phot', '--noned', '--noplot']

    path="tests/"

    filephot = "imgblock-out.txt"
    filephot= path+filephot

    filecomp = "imgblock-comp.fits"
    filecomp= path+filecomp

    filepng = "imgblock.png"
    filemulpng = "imgblock-mul.png"

    filesig = "imgblock-sig.fits"
    filesig= path+filesig

    filecheck = "imgblock-check.fits"
    filecheck= path+filecheck

    filecube = "imgblock-cub.png"
    filecube = path+filecube



    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)




    # read user's input 
    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)




    #tolerance parameter
    tol = 1e-3

    bt= 1.000 
    tidal = 2.011 
    lchinu = 1.043 
    bump = 0.129 
    snr = 1.477 
    std_snr = 3.860 
    aic= 1399.730 
    bic = 1436.107 

    diffbt =abs(bt-photapi.BulgeToTotal ) 
    difftidal = abs(tidal-photapi.tidal)
    difflchinu = abs(lchinu-photapi.objchinu)
    diffbump = abs(bump-photapi.bump)
    diffsnr = abs(snr-photapi.snr)
    diffstd_snr = abs(std_snr-photapi.stdsnr)
    diffaic= abs(aic-photapi.AICrit)
    diffbic = abs(bic-photapi.BICrit)



    assert os.path.isfile(filephot)

    assert diffbt < tol
    assert difftidal < tol
    assert difflchinu < tol
    assert diffbump < tol
    assert diffsnr < tol
    assert diffstd_snr < tol
    assert diffaic < tol
    assert diffbic < tol



    if os.path.isfile(filepng):
        os.remove(filepng)

    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)


# checking the computation of sky with the gradient method
def test_gsky():


    arg=['tests/galfit.01','-gsky', '-ri','2', '-skw','3','--noplot']

    path="tests/"

    filephot = "imgblock-out.txt"
    filephot= path+filephot

    filecomp = "imgblock-comp.fits"
    filecomp= path+filecomp

    filepng = "imgblock.png"
    filemulpng = "imgblock-mul.png"

    filesig = "imgblock-sig.fits"
    filesig= path+filesig

    filecheck = "imgblock-check.fits"
    filecheck= path+filecheck

    filecube = "imgblock-cub.png"
    filecube = path+filecube

    filersky = "imgblock-ring.fits"
    filersky = path + filersky

    filermsky = "imgblock-ringmask.fits"
    filermsky = path + filermsky



    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)
 


    # read user's input 
    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)




    #tolerance parameter
    tol = 1e-2

    sky = 0.63


    diffsky = abs(sky - photapi.gradskymean)



    assert diffsky < tol



    if os.path.isfile(filepng):
        os.remove(filepng)

    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)

 # checking the computation of sky with the random method
def test_rsky():


    arg=['tests/galfit.01','-rsky','--noplot']

    path="tests/"

    filephot = "imgblock-out.txt"
    filephot= path+filephot

    filecomp = "imgblock-comp.fits"
    filecomp= path+filecomp

    filepng = "imgblock.png"
    filemulpng = "imgblock-mul.png"

    filesig = "imgblock-sig.fits"
    filesig= path+filesig

    filecheck = "imgblock-check.fits"
    filecheck= path+filecheck

    filecube = "imgblock-cub.png"
    filecube = path+filecube

    filersky = "imgblock-ring.fits"
    filersky = path + filersky

    filermsky = "imgblock-ringmask.fits"
    filermsky = path + filermsky



    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)
 


    # read user's input 
    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)




    #tolerance parameter
    tol = 5

    sky = 0.9


    diffsky = abs(sky - photapi.randskymean)



    assert diffsky < tol



    if os.path.isfile(filepng):
        os.remove(filepng)

    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)
 
