cones
=====

Scattering cone brightness fitting tools

I have made this repository to hold the code that I used to fit the lateral and radial scattering cone profiles. I will update it regularly with new pieces of code and modules in the hope that they may be reused later.

The following is a description of the contents of each .py file.


constants.py
=============================
This file holds physical constants in appropriate units for use in the other files

smc_dp.py
============================
This is a module that reads and extrapolates the Draine phase functions for use in obtaining the differential cross-section for dust. This is essential in the fitting of a dust-only model.

smc_pol.py
============================
This is a module that reads the numerical polarization fraction functions and extrapolates to obtain piecewise linear functions over the entire range of angles.


calc.py
===========================
This is an executabel file that calculates the properties of multiple objects and outputs a summary of these properties as well as normalized radial and lateral sections. The summary includes luminosity (L), distances, etc. The input to this file is a list of files obtained from reading brightness profiles in ds9. These files also include an additional header that specifies the length (angular dimension) of the sections taken and their position.

sim_ele_psf.py
==========================
This is the main executable file that calculates the ultimate functions to fit the profiles taking into account the PSF. The input to this executable is a file (called inits.txt) holding the initial values of the fit parameters as well as a collection of files containing the normalized lateral and radial profiles. The normalized profiles are produced as the output of calc.py when run on 'raw' sections like the example provided. It takes several options when being run and these are illustrated in the following examples:

$ python sim_ele_psf.py --plotInits all
This command reads a file called inits.txt and uses the values of the parameters to make plots of all objects. No fitting procedure is performed

$ python sim_ele_psf.py --display --fit 2,3,4
This command reads the initial value of the parameters and fits the radial and lateral profiles independently. This is done for objects 2, 3 and 4 only. The --display option instructs the program to show the plots as they are being made and asks for user input to continue.

$ python sim_ele_psf.py --sim --fit 10
This command performs a simultaneous fit of the lateral and radial profiles again using initial values from the file inits.txt.

All these fits are electron-only models.

