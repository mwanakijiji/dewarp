# This takes a standalone list of polynomial coefficients and makes a dewarp plot
# to visualize the distortion, in much the same way as the script find_dewarp_solution.py
# does. (Note, however, that some of the 'pre-' and 'post-' conventions are flipped
# here, due to the way the plot is made in find_dewarp_solution.py, and the nontriviality
# of finding inverses to Kx and Ky (which are essentially multidimensional matrices) from each other.) 

# created 23 Dec 2016 by E.S.

import numpy as np
from astrom_lmircam_soln import *
from astrom_lmircam_soln import dewarp
from astrom_lmircam_soln import make_barb_plot
import matplotlib.pyplot as plt



#####################################################################
# FIND THE MAPPING BETWEEN THE IDEAL AND EMPIRICAL PINHOLE COORDS 

# coefficients from the LEECH pipeline based on Maire+ 2015 measurements
Kx = [[-2.1478925,    0.0058138110,  -7.6396687e-06,   2.5359596e-09],
      [1.0109149,  -2.3826537e-05,   2.8458629e-08,  -9.3206482e-12],
      [-2.1164521e-05,   5.3115381e-08,  -6.6315643e-11,   2.2888432e-14],
      [1.2983972e-08,  -4.1253977e-11,   5.1637044e-14,  -1.5988376e-17]]
Ky = [[9.2717864,      0.98776733,   4.3514612e-06,   9.3450739e-09],
      [-0.013617797,  -3.9526096e-05,   8.1204222e-08,  -5.2048768e-11],
      [1.1313247e-05,   6.7127301e-08,  -1.6531988e-10,   1.0656544e-13],
      [1.6283111e-09,  -2.7723216e-11,   8.2118035e-14,  -5.3695050e-17]]


# make an even grid to tack onto the POST-dewarped image; these grid points are the locations where the barbs will be
x_vec_grid_postdewarp, y_vec_grid_postdewarp, xIntervalsOnly, yIntervalsOnly = make_barb_plot.barb_grid()


# map the points from the post-dewarp barb grid to the pre-dewarp grid 
predewarp_coords = dewarp.make_dewarp_coordinates_sparseGrid(xIntervalsOnly,
                                                             yIntervalsOnly,
                                                             np.array(Kx).T,
                                                             np.array(Ky).T)

# remove stray array dimensions of thickness 1
x_vec_grid_predewarp = np.squeeze(predewarp_coords)[1,:,:]
y_vec_grid_predewarp = np.squeeze(predewarp_coords)[0,:,:]


# find the difference between the pre- and post-dewarp grids
x_vec_grid_diff = np.subtract(x_vec_grid_postdewarp,x_vec_grid_predewarp)
y_vec_grid_diff = np.subtract(y_vec_grid_postdewarp,y_vec_grid_predewarp)


# set up the barb (quiver) grids
barbs = plt.quiver(np.ravel(x_vec_grid_predewarp), np.ravel(y_vec_grid_predewarp), # barbs in the sampled region
                   np.ravel(x_vec_grid_diff), np.ravel(y_vec_grid_diff),color='k',scale=0.3, scale_units='xy')


# finally, make the barb plot
plt.axis([0.,2055.,0.,2055.])
ax = plt.axes()
ax.set_aspect(1.)
axes = plt.gca()
axes.set_xlim([0,2050])
plt.title('LMIRCam Distortion Correction Map\n(epoch between July and Dec 2014; see Maire+ 2015 A&A 576 A133)')
plt.xlabel('x (pix)')
plt.ylabel('y (pix)')
plt.tight_layout()
plt.show()
