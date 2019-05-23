###########################

# This is the script to run for making a warping/astrometric solution to LMIRCam, using data taken
# in Nov and Dev 2016

# parent find_dewarp_soln.py created by E.S., Nov 2016
# child apply_dewarp_soln.py made by E.S., Feb 2017
# edits by E.S., May 2019

import numpy as np
from astrom_lmircam_soln import *
from astrom_lmircam_soln import dewarp
from astropy.io import fits


#####################################################################
# SET THE DEWARP COEFFICIENTS

# dewarp coefficients

# SX, UT 2019 Jan 24:
'''
Kx = [[-1.54869268e+01,  2.17731891e-02, -1.16311438e-05,  1.72268320e-09],
      [ 1.03589421e+00, -2.52297644e-05,  2.25845524e-08, -6.82857090e-12],
      [-2.73357199e-05,  8.26938031e-09, -1.50590654e-11,  7.29255070e-15],
      [ 8.40890698e-09, -2.20628759e-12,  4.65179926e-15, -2.36037709e-18]]

Ky = [[ 6.81966212e+00,  9.95703945e-01, -7.47825973e-06,  8.28967040e-09],
      [-1.98805359e-02, -2.80102218e-05,  9.51421116e-09, -2.96630491e-12],
      [ 5.78674771e-06,  2.15682679e-08, -1.35997721e-11,  4.75076932e-15],
      [ 1.60545966e-09, -5.27273081e-12,  5.72485139e-15, -2.06441615e-18]]
'''
      
# DX, UT 2019 Jan 24:
'''
Kx = [[-1.48945923e+01,  2.00898026e-02, -9.60604219e-06, 9.29783279e-10],
      [ 1.03374217e+00, -1.86557902e-05,  1.37587067e-08, -3.08994609e-12],
      [-2.50236726e-05,  8.91493383e-10, -4.44492395e-12, 2.66171048e-15],
      [ 7.63440010e-09,  3.87362904e-13,  8.22451038e-16, -6.85924248e-19]]

Ky = [[7.29742888e+00,  9.94071005e-01, -5.83084023e-06, 7.76900384e-09],
      [-2.24096548e-02, -1.93099545e-05,  3.87183244e-10, 9.34686744e-14],
      [ 8.92216788e-06,  1.03140645e-08, -1.40953784e-12, 5.64959139e-16],
      [ 5.44087872e-10, -1.29336418e-12,  1.27236764e-15, -5.01793976e-19]]
'''
      
#####################################################################
# DEWARP TEST FILE

# this is just to get the shape; image content doesn't matter
hdul = fits.open('test_pre_dewarping.fits')
imagePinholes = hdul[0].data.copy()

# map the coordinates that define the entire image plane
# (transposition due to a coefficient definition change btwn Python and IDL)
dewarp_coords = dewarp.make_dewarp_coordinates(imagePinholes.shape,
                                               np.array(Kx).T,
                                               np.array(Ky).T) 


## CAN START FOR-LOOP HERE

# grab the pre-dewarp image and header
imageAsterism, header = fits.getdata('test_pre_dewarping.fits',
                                     0, header=True)

# dewarp the image
dewarpedAsterism = dewarp.dewarp_with_precomputed_coords(imageAsterism,
                                                         dewarp_coords,
                                                         order=3)

# write out
fits.writeto('test_post_dewarping.fits',
             np.squeeze(dewarpedAsterism),
             header,
             overwrite=False)

#####################################################################
