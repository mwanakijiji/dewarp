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
from pytictoc import TicToc
import regions
from regions.core import PixCoord
from regions.shapes.circle import CirclePixelRegion
import configparser
import glob

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
config.read("astrom/config.ini")

#####################################################################
# FIND THE MAPPING BETWEEN THE IDEAL AND EMPIRICAL PINHOLE COORDS 

# retrieve pickled pinhole coordinates ## this fcn def is repeated from previous module
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


# retrieve pickled pinhole coordinates
def get_dewarp_soln(fileString='test'):
    '''
    read the matched and sorted pickled pinhole coordinates
    '''
    picklefile = 'dewarp_soln_'+fileString+'.pkl'
    fo1=open(picklefile,'rb')
    dat=pickle.load(fo1)
    fo1.close()
    Kx = np.array(dat['Kx']) # FYI
    Ky = np.array(dat['Ky']) # FYI
    dewarp_coords =  np.array(dat['dewarp_coords']) # this is what the pipeline actually uses
    return dewarp_coords

def apply_dewarp(writeoutString='test',maskUnsampled=False):
    
    #####################################################################
    # APPLY THE DEWARP SOLUTION TO THE SCIENCE IMAGES

    # obtain the list of files to dewarp
    asterism_frames_directory_retrieve = str(config["data_dirs"]["DIR_ASTERISM_BASIC"])
    asterism_frames_pre_dewarp_names = list(glob.glob(os.path.join(asterism_frames_directory_retrieve, "*.fits")))

    # retrieve dewarp solution
    dewarp_coords = get_dewarp_soln(writeoutString)

    # option: mask region unsampled by pinholes
    if maskUnsampled:
        print('------------------------------')
        print('Calculating mask for unsampled region...')
        xy_model_d, xy_model_not_d, xy_empirical = get_pinhole_grid_coords(writeoutString) # retrieve pinhole locations

        pixcoord = PixCoord(xy_empirical[:,0],xy_empirical[:,1]) # the centers of regions around which to draw radii
        radialDist = 100 # consider pixels within this radius of any pinhole to be 'sampled'
        regions = [CirclePixelRegion(center=PixCoord(x,y), radius=radialDist) for x,y in xy_empirical] # overlap the circles centered on each pinhole

        ## need to read in array shape in argument below
        mask_canvas = np.ones((2048,2048), dtype=np.int32) # initialize mask of ones
            
        for x,y in xy_empirical:
            region_i = CirclePixelRegion(center=PixCoord(x,y), radius=radialDist)
            mask_i = region_i.to_mask()
            ## need to read in array shape in argument below
            thisPatch = np.add(1.,-mask_i.to_image((2048,2048))).astype(np.int32) # apply mask to a larger canvas
            mask_canvas = np.multiply(mask_canvas,thisPatch).astype(np.int32) # add patch

        mask_final = np.add(1.,-mask_canvas).astype(np.int32) # 1=sampled, 0=unsampled
    
    print('------------------------------')
    for frameNum in range(0,len(asterism_frames_pre_dewarp_names)):
        print('Dewarping frame %05i'%frameNum+'...')

        t = TicToc() # create instance of timer
        t.tic() # start timer

        # grab the pre-dewarp image and header
        imageAsterism, header = fits.getdata(asterism_frames_pre_dewarp_names[frameNum],
                                             0, header=True)

        # option: apply mask
        if maskUnsampled:
            imageAsterism = np.add(
                np.multiply(np.squeeze(imageAsterism),mask_final),
                mask_canvas*np.median(imageAsterism[550:560,550:560]
                )).astype(np.int32) # mask the unsampled region with a median value
        
        # dewarp the image
        dewarpedAsterism = dewarp.dewarp_with_precomputed_coords(imageAsterism,
                                                             dewarp_coords,
                                                             order=3)
        
                    
        # write out
        fits.writeto(str(config["data_dirs"]["DIR_ASTERISM_DEWARP"] + \
                         "dewarped_" + \
                         os.path.basename(asterism_frames_pre_dewarp_names[frameNum])),
                     np.squeeze(dewarpedAsterism), header, overwrite=False)

        t.toc()
        print('------------------------------')
