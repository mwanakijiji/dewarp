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
from scipy import spatial
import configparser

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
config.read("astrom/config.ini")


#####################################################################
# FIND THE MAPPING BETWEEN THE IDEAL AND EMPIRICAL PINHOLE COORDS 

# retrieve pickled pinhole coordinates
def get_pinhole_grid_coords(fileString):
    '''
    read the matched and sorted pickled pinhole coordinates
    '''
    picklefile = 'pinhole_coordinates_'+fileString+'.pkl'
    fo1=open(picklefile,'rb')
    dat=pickle.load(fo1)
    fo1.close()
    coordsModel_d = np.array(dat['coordsModel_d_matched']) # not sure what this would be useful for, though
    coordsModel_not_d = np.array(dat['coordsModel_not_d_matched']) # 'ideal' (x,y)
    coordsEmpirical =  np.array(dat['coordsEmpirical_matched']) # empirical (x,y)
    return coordsModel_d, coordsModel_not_d, coordsEmpirical

def find_dewarp(fileString='test',dateString='---',plot=True):
    ## need to move this elsewhere

    # map the x,y lists together using a fcn based directly off the IDL fcn polywarp.pro (I tested it against another Python version (polywarp_v2.py), too)
    # note order of the arguments: adopting the IDL documentation notation on polywarp.pro, the warped coords are (xi,yi) and dewarped coords (xo,yo)
    xy_model_d, xy_model_not_d, xy_empirical = get_pinhole_grid_coords(fileString=fileString)

    # get coefficients
    Kx, Ky = polywarp.polywarp(xy_empirical,
                           xy_model_not_d,
                           degree=3)


    # retrieve array for checking its size
    hdul = fits.open(config["data_dirs"]["DIR_PINHOLE_BASIC"] + config["src_file_names"]["PINHOLE_FITS"]) 
    imagePinholes = hdul[0].data.copy()

    # map the coordinates that define the entire image plane
    # note the below couple functions appear in the LEECH pipeline, and the above Kx, Ky are oriented in the same convention as the Kx, Ky in the LEECH pipeline
    dewarp_coords = dewarp.make_dewarp_coordinates(imagePinholes.shape,
                                               np.array(Kx).T,
                                               np.array(Ky).T) # transposed due to a coefficient definition change btwn Python and IDL

    # print to screen
    print('Dewarp coefficients:')
    print('Kx:')
    print(Kx)
    print('Ky:')
    print(Ky)
    print('------------------------------')
    print('Resampling the array...')
    print('------------------------------')

    # pickle the solution
    picklename='dewarp_soln_'+fileString+'.pkl'
    fo=open(picklename,'wb')
    pickle.dump({'Kx':Kx,
             'Ky':Ky,
             'dewarp_coords':dewarp_coords},
            fo)
    fo.close()
        
    print('Saved dewarp solution to:')
    print(picklename)
    print('------------------------------')


    #####################################################################
    # MAKE A DEWARP VECTOR PLOT TO VISUALIZE THE DISTORTION
    # vectors indicate 'this is the direction a point here should be shifted'

    if plot:
        # make an even grid to tack onto the pre-dewarped image; these grid points are the locations where the barbs will be
        x_vec_grid_predewarp, y_vec_grid_predewarp, xIntervalsOnly, yIntervalsOnly = make_barb_plot.barb_grid()

        # do the inverse fit, to find the coefficients that take us from warped coords --> dewarped coords
        KxInv, KyInv = polywarp.polywarp(xy_model_not_d,
                                 xy_empirical,
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
        # this can be tuned as you see fit
        N = 36. 

        # ... find which barbs which are in the 'sampled' region and which are in the 'projected' region

        closeOnes, farOnes = make_barb_plot.find_close_pts_cdist(np.transpose([np.ravel(x_vec_grid_predewarp),
                                                                       np.ravel(y_vec_grid_predewarp)]),
                                                         np.transpose([xy_empirical[:,0],
                                                                       xy_empirical[:,1]]),
                                                         N)

        # set up the barb (quiver) grids
        sampledReg = plt.quiver((np.ravel(x_vec_grid_predewarp))[closeOnes], (np.ravel(y_vec_grid_predewarp))[closeOnes], # barbs in the sampled region
           (np.ravel(x_vec_grid_diff))[closeOnes], (np.ravel(y_vec_grid_diff))[closeOnes],color='k')
        unsampledReg = plt.quiver((np.ravel(x_vec_grid_predewarp))[farOnes], (np.ravel(y_vec_grid_predewarp))[farOnes], # barbs in the projected region
(np.ravel(x_vec_grid_diff))[farOnes], (np.ravel(y_vec_grid_diff))[farOnes],color='grey',width=0.001)

        # add key
        sampledKey = plt.quiverkey(sampledReg, 1.2, 0.8, 12, 'sampled region', labelpos='N',\
                           fontproperties={'weight': 'normal','size': 'medium'})
        unsampledKey = plt.quiverkey(unsampledReg, 1.2, 0.6, 12, r'projection in' '\n' r'unsampled region', labelpos='N',\
                           fontproperties={'weight': 'normal','size': 'medium'})

        # finally, make the barb plot
        plt.axis([0.,2055.,0.,2055.])
        ax = plt.axes()
        ax.set_aspect(1.)
        axes = plt.gca()
        axes.set_xlim([0,2050])
        plt.title('LMIRCam Distortion Correction Map\n(based on pinhole data taken '+dateString+')')
        plt.xlabel('x (pix)')
        plt.ylabel('y (pix)')
        plt.tight_layout()
        print('Close the plot to proceed...')
        plt.savefig("barb_plot.pdf", overwrite=True)
        plt.show()
