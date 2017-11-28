# This finds star locations in pixel space in asterism images that have been dewarped (using the
# solution found by find_dewarp_solution.py) and derotated (based on the PA in the FITS header).
# The coordinates (including false positives) are printed to the Terminal and have to be manually
# cross-checked with known star positions, after which the pixel coordinates have to be written
# into the appropriate dictionaries in find_plate_scale_and_orientation.py.

# I know it's inelegant and a lot of manual work-- if you know a better way, implement it!

# created Dec 2016 by E.S.

import numpy as np
from astrom_lmircam_soln import *
from astrom_lmircam_soln import find_pinhole_centroids
import scipy
from astropy.io import fits
import matplotlib.pyplot as plt
import photutils
from photutils import DAOStarFinder

number_of_dithers = 18 # number of dither positions over the asterism
star_coords_every_dither = {} # initialize dictionary to hold the coords of found stars

for ditherPos in range(0,number_of_dithers):
    print("Finding star positions in dither position "+str(ditherPos)+" ...")

    # read in image and header
    imageMedian, header = fits.getdata(calibrated_trapezium_data_stem+
                                       'step04_ditherMedians/median_dither_'+'%02i'%ditherPos+'.fits',
                                       0, header=True)
    
    keyName = "dither_pos_"+str(ditherPos)

    # find star locations; input parameters may need some find-tuning to get good results
    xCoordsAsterism, yCoordsAsterism = find_pinhole_centroids.find_psf_centers(imageMedian,20.,15000.) 
    star_coords_every_dither[keyName] = np.transpose([xCoordsAsterism, yCoordsAsterism])
    print(np.transpose([xCoordsAsterism,yCoordsAsterism]))
    print("Please check the plot and note true positives among the printed star positions. Close the plot to continue.")
    
    # plots for cross-checking found centroids with true stars (save plots manually as desired)
    plt.imshow(imageMedian, origin="lower")
    plt.scatter(xCoordsAsterism, yCoordsAsterism)
    plt.title("LMIRCam Trapezium observation, UT 2016 Nov 12\ndither position "+"%02i"%ditherPos)
    plt.show()
    
# optional: save centroids as a FYI pickle binary file (remember: some coordinates are false positives! need to manually check)
#pickle.dump(star_coords_every_dither, open("centroid_coords_trapezium_ut_2016_nov_12.p", "wb"))
