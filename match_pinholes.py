# This is the script to run for making a warping/astrometric solution to LMIRCam, using data taken
# in Nov and Dev 2016

# created by E.S., Nov 2016
# revamped to be more user-friendly, Nov 2017

import numpy as np
from astrom_lmircam_soln import *
from astrom_lmircam_soln import find_pinhole_centroids
from astrom_lmircam_soln import polywarp
from astrom_lmircam_soln import polywarp_v2 # the check for polywarp
from astrom_lmircam_soln import dewarp
from astrom_lmircam_soln import make_barb_plot
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import ipdb


#####################################################################
# RETRIEVE MEDIAN PINHOLE GRID IMAGE, FIND PINHOLE CENTERS

# write the git hash
write_hash('git_hash_match_pinholes.txt')

# (note this section will require multiple run-throughs until optimal parameters are specified: the ideal grid coords, the missed pinhole coords, etc.)

hdul = fits.open(calibrated_pinholes_data_stem+'step02_medianed/pinhole_image_median_DXonly_171108.fits') # median image of pinholes, with vignetted regions masked
imagePinholes = hdul[0].data.copy()

xCoordsIdealFullGrid, yCoordsIdealFullGrid = find_pinhole_centroids.put_down_grid_guesses([100,-100],48.0,[512,512],2e-8,0.65) # sets down an 'ideal' set of pinholes made to match the real pinholes as closely as possible
# approxHoleSpacingPass,barrelCenterPass,barrelDegreePass,rotationAnglePass
import ipdb; ipdb.set_trace()
xCoordsFoundAutomated, yCoordsFoundAutomated = find_pinhole_centroids.find_psf_centers(imagePinholes,20.,50000.) # finds the actual st of pinholes

# manually-found locations of pinholes that the above routine missed
xCoordsMissed = [71.774,697.353,1460.66]
yCoordsMissed = [1267.57,1404.06,737.932]

# remove subset of pinhole guesses that are in the vignetted region
xCoordsIdeal = xCoordsIdealFullGrid[np.where(np.logical_or(xCoordsIdealFullGrid>36,
                                                               yCoordsIdealFullGrid>760))]
yCoordsIdeal = yCoordsIdealFullGrid[np.where(np.logical_or(xCoordsIdealFullGrid>36,
                                                               yCoordsIdealFullGrid>760))]
ipdb.set_trace()
# concatenate arrays
xCoordsFound = np.concatenate((xCoordsFoundAutomated,
                               xCoordsMissed),
                              axis=0)
yCoordsFound = np.concatenate((yCoordsFoundAutomated,
                               yCoordsMissed),
                              axis=0)

# sort the arrays to make each set of ideal and empirical x,y entries correspond to the same pinhole
xIdeal_sorted, yIdeal_sorted, xFound_sorted, yFound_sorted = find_pinhole_centroids.consistency_xy_found_guesses(xCoordsIdeal,
                                                                                                                 yCoordsIdeal,
                                                                                                                 xCoordsFound,
                                                                                                                 yCoordsFound)
