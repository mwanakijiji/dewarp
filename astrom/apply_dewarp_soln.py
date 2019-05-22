# This is the script to run for making a warping/astrometric solution to LMIRCam, using data taken
# in Nov and Dev 2016

# parent find_dewarp_soln.py created by E.S., Nov 2016
# child apply_dewarp_soln.py made by E.S., Feb 2017

import numpy as np
from astrom_lmircam_soln import *
from astrom_lmircam_soln import polywarp
from astrom_lmircam_soln import dewarp
from astropy.io import fits
import matplotlib.pyplot as plt


#####################################################################
# SET THE DEWARP COEFFICIENTS

# the below coefficients are relevant for LMIRCam 2048x2048 readouts after modifications in summer 2016
# (in Maire+ 2015 format, here x_0 = y_0 = 0)
Kx = [[ -4.74621436e+00,   9.99369200e-03,  -4.69741638e-06,   4.11937105e-11],
      [  1.01486148e+00,  -2.84066638e-05,   2.10787962e-08,  -3.90558311e-12],
      [ -1.61139243e-05,   2.24876212e-08,  -2.29864156e-11,   6.59792237e-15],
      [  8.88888428e-09,  -1.03720381e-11,   1.05406782e-14,  -3.06854175e-18]]

Ky = [[  9.19683947e+00,   9.84613002e-01,  -1.28813904e-06,   6.26844974e-09],
      [ -7.28218373e-03,  -1.06359740e-05,   2.43203662e-09,  -1.17977589e-12],
      [  9.48872623e-06,   1.03410741e-08,  -2.38036199e-12,   1.17914143e-15],
      [  3.56510910e-10,  -5.62885797e-13,  -5.67614656e-17,  -4.21794191e-20]]

#####################################################################
# DEWARP TEST FILE

hdul = fits.open('pinhole_image_median_vignettingBlocked.fits') # this is just to get the shape
imagePinholes = hdul[0].data.copy()

# map the coordinates that define the entire image plane
dewarp_coords = dewarp.make_dewarp_coordinates(imagePinholes.shape,
                                               np.array(Kx).T,
                                               np.array(Ky).T) # transposed due to a coefficient definition change btwn Python and IDL


# grab the pre-dewarp image and header
imageAsterism, header = fits.getdata('pinhole_image_median_vignettingBlocked.fits',
                                     0, header=True)

# dewarp the image
dewarpedAsterism = dewarp.dewarp_with_precomputed_coords(imageAsterism,
                                                         dewarp_coords,
                                                         order=3)

# write out
fits.writeto('test_dewarped.fits',
             np.squeeze(dewarpedAsterism),
             header,
             clobber=False)

#####################################################################
