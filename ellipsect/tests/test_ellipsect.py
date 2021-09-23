

import pytest
import os


from ellipsect import ellipsectors 

from ellipsect.inout.read import InputSys
from ellipsect.sectors.sect import SectorsGalfit


# simple run test
def test_exit():
    with pytest.raises(SystemExit) as e:
        ellipsectors.run()
    assert e.type == SystemExit 
    assert e.value.code == 0 

# checking the creation of files
def test_files():

    arg=['./ellsec.py', 'ellipsect/tests/galfit.01', '-noplot']

    path="ellipsect/tests/"

    filepng = "imgblock.png"
    filemulpng = "imgblock-mul.png"

    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)


    params = InputSys(arg)
    photapi = SectorsGalfit(params)

    assert os.path.isfile(filepng)
    assert os.path.isfile(filemulpng)


    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)



  

# checking the creation of the components file
def test_comp():

    arg=['./ellsec.py', 'ellipsect/tests/galfit.01','-comp', '-noplot']

    path="ellipsect/tests/"

    filepng = "imgblock.png"
    filemulpng = "imgblock-mul.png"

    filepng = path+filepng
    filemulpng = path+filemulpng



    filecomp = "imgblock-comp.fits"

    filecomp= path+filecomp

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    params = InputSys(arg)
    photapi = SectorsGalfit(params)

    assert os.path.isfile(filecomp),"is GALFIT installed?"


    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)


    if os.path.isfile(filecomp):
        os.remove(filecomp)







# checking the creation of the sbout files
def test_phot():


    arg=['./ellsec.py', 'ellipsect/tests/galfit.01','-phot', '-noned', '-noplot']

    path="ellipsect/tests/"

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




    params = InputSys(arg)
    photapi = SectorsGalfit(params)

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



