# This is the script to run for making a warping/astrometric solution to LMIRCam, using data taken
# in Nov and Dev 2016

# created by E.S., Nov 2016

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


#####################################################################
# RETRIEVE MEDIAN PINHOLE GRID IMAGE, FIND PINHOLE CENTERS

# (note this section will require multiple run-throughs until optimal parameters are specified: the ideal grid coords, the missed pinhole coords, etc.)

hdul = fits.open(calibrated_pinholes_data_stem+'pinhole_image_median_vignettingBlocked.fits') # median image of pinholes, with vignetted regions masked
imagePinholes = hdul[0].data.copy()
xCoordsIdealFullGrid, yCoordsIdealFullGrid = find_pinhole_centroids.put_down_grid_guesses(48.0,0.65) # sets down an 'ideal' set of pinholes made to match the real pinholes as closely as possible
xCoordsFoundAutomated, yCoordsFoundAutomated = find_pinhole_centroids.find_psf_centers(imagePinholes,20.,50000.) # finds the actual st of pinholes

# manually-found locations of pinholes that the above routine missed
xCoordsMissed = [71.774,697.353,1460.66]
yCoordsMissed = [1267.57,1404.06,737.932]

# remove subset of pinhole guesses that are in the vignetted region
xCoordsIdeal = xCoordsIdealFullGrid[np.where(np.logical_or(xCoordsIdealFullGrid>36,
                                                               yCoordsIdealFullGrid>760))]
yCoordsIdeal = yCoordsIdealFullGrid[np.where(np.logical_or(xCoordsIdealFullGrid>36,
                                                               yCoordsIdealFullGrid>760))]

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

#####################################################################
# FIND THE MAPPING BETWEEN THE IDEAL AND EMPIRICAL PINHOLE COORDS 

# map the x,y lists together using a fcn based directly off the IDL fcn polywarp.pro (I tested it against another Python version (polywarp_v2.py), too)
# note order of the arguments: adopting the IDL documentation notation on polywarp.pro, the warped coords are (xi,yi) and dewarped coords (xo,yo)
Kx, Ky = polywarp.polywarp(xFound_sorted,
                           yFound_sorted,
                           xIdeal_sorted,
                           yIdeal_sorted,
                           degree=3)

# map the coordinates that define the entire image plane
# note the below couple functions appear in the LEECH pipeline, and the above Kx, Ky are the same as the Kx, Ky in the LEECH pipeline
dewarp_coords = dewarp.make_dewarp_coordinates(imagePinholes.shape,
                                               np.array(Kx).T,
                                               np.array(Ky).T) # transposed due to a coefficient definition change btwn Python and IDL

# optional: view the dewarped pinhole image as a check
#dewarpedImg = dewarp.dewarp_with_precomputed_coords(imagePinholes,dewarp_coords,order=3) # np.squeeze shouldn't be a problem for displaying, right?
#plt.imshow(dewarpedImg)
#plt.show()


#####################################################################
# APPLY THE DEWARP SOLUTION TO THE SCIENCE IMAGES 
'''
for frameNum in range(1892,2252):
    print('Dewarping frame %05i'%frameNum+'...')

    # grab the pre-dewarp image and header
    imageAsterism, header = fits.getdata(calibrated_trapezium_data_stem+
                                       'step01_darkSubtBadPixCorrect/lm_161112_'+'%05i'%frameNum+'.fits',
                                       0, header=True)

    # dewarp the image
    dewarpedAsterism = dewarp.dewarp_with_precomputed_coords(imageAsterism,
                                                             dewarp_coords,
                                                             order=3)

    # write out
    fits.writeto(calibrated_trapezium_data_stem+
                 'step02_dewarped/lm_161112_'+'%05i'%frameNum+'.fits',
                 np.squeeze(dewarpedAsterism), header, clobber=True)
'''

#####################################################################
# MAKE A DEWARP VECTOR PLOT TO VISUALIZE THE DISTORTION


# make an even grid to tack onto the pre-dewarped image; these grid points are the locations where the barbs will be
x_vec_grid_predewarp, y_vec_grid_predewarp, xIntervalsOnly, yIntervalsOnly = make_barb_plot.barb_grid()

# do the inverse fit, to find the coefficients that take us from warped coords --> dewarped coords
KxInv, KyInv = polywarp.polywarp(xIdeal_sorted,
                                 yIdeal_sorted,
                                 xFound_sorted,
                                 yFound_sorted,
                                 degree=3)

# map the points from the pre-dewarp barb grid to the post-dewarp grid 
dewarp_coords_inv = dewarp.make_dewarp_coordinates_sparseGrid(xIntervalsOnly,
                                                              yIntervalsOnly,
                                                              np.array(KxInv).T,
                                                              np.array(KyInv).T)

# remove stray array dimensions of thickness 1
x_vec_grid_postdewarp = np.squeeze(dewarp_coords_inv)[1,:,:]
y_vec_grid_postdewarp = np.squeeze(dewarp_coords_inv)[0,:,:]

# find the difference between the pre- and post-dewarp grids
x_vec_grid_diff = np.subtract(x_vec_grid_postdewarp,x_vec_grid_predewarp)
y_vec_grid_diff = np.subtract(y_vec_grid_postdewarp,y_vec_grid_predewarp)

# set cutoff distance between barb points and the pinhole locations that were actually sampled;
# this is needed to ...
N = 25. 

# ... find which barbs which are in the 'sampled' region and which are in the 'projected' region
closeOnes, farOnes = make_barb_plot.find_close_pts_cdist(np.transpose([np.ravel(x_vec_grid_predewarp),
                                                                       np.ravel(y_vec_grid_predewarp)]),
                                                         np.transpose([xFound_sorted,
                                                                       yFound_sorted]),
                                                         N)

# set up the barb (quiver) grids
sampledReg = plt.quiver((np.ravel(x_vec_grid_predewarp))[closeOnes], (np.ravel(y_vec_grid_predewarp))[closeOnes], # barbs in the sampled region
                        (np.ravel(x_vec_grid_diff))[closeOnes], (np.ravel(y_vec_grid_diff))[closeOnes],color='k') # ,scale=0.3, scale_units='xy') # can fix the scale of the arrow lengths
unsampledReg = plt.quiver((np.ravel(x_vec_grid_predewarp))[farOnes], (np.ravel(y_vec_grid_predewarp))[farOnes], # barbs in the projected region
                          (np.ravel(x_vec_grid_diff))[farOnes], (np.ravel(y_vec_grid_diff))[farOnes],color='grey') # ,scale=0.3, scale_units='xy')

# add key
sampledKey = plt.quiverkey(sampledReg, 1.2, 0.8, 12, 'sampled region', labelpos='N',\
                           fontproperties={'weight': 'normal','size': 'medium'})
unsampledKey = plt.quiverkey(unsampledReg, 1.2, 0.6, 12, r'projection in' '\n' r'unsampled region' '\n' r'(mostly vignetted)', labelpos='N',\
                           fontproperties={'weight': 'normal','size': 'medium'})

# finally, make the barb plot
plt.axis([0.,2055.,0.,2055.])
ax = plt.axes()
ax.set_aspect(1.)
# optional: plate scale and orientation annotations, which have to be determined with find_plate_scale_and_orientation.py, using dewarped asterism images
''' 
plateScale_ann = ax.annotate("plate scale 10.69 $\pm$0.04 mas/pix\n(based on Trapezium data\ntaken UT 2016 Nov 12)",
                  xy=(1.2, 0.4), xycoords='axes fraction',
                  xytext=(1.2, 0.4), textcoords='axes fraction',
                  size=12, va="center", ha="center",
                  )
orientation_ann = ax.annotate("for true North orientation,\nrotate dewarped, derotated\nimages by $1.1^{\circ}$ W of N\n(also based on Trapezium data\ntaken UT 2016 Nov 12)",
                  xy=(1.2, 0.4), xycoords='axes fraction',
                  xytext=(1.2, 0.2), textcoords='axes fraction',
                  size=12, va="center", ha="center",
                  )
'''
axes = plt.gca()
axes.set_xlim([0,2050])
plt.title('LMIRCam Distortion Correction Map\n(based on pinhole data taken 2016 Nov 22)')
plt.xlabel('x (pix)')
plt.ylabel('y (pix)')
plt.tight_layout()
plt.show()
