# This finds star locations in pixel space in asterism images that have been dewarped (using the
# solution found by find_dewarp_solution.py) and derotated (based on the PA in the FITS header).
# The coordinates (including false positives) are printed to the Terminal and have to be manually
# cross-checked with known star positions, after which the pixel coordinates have to be written
# into the appropriate dictionaries in find_plate_scale_and_orientation.py.

# created Dec 2016 by E.S.

import numpy as np
from astrom_lmircam_soln import *
from astrom_lmircam_soln import find_pinhole_centroids
import scipy
from astropy.io import fits
import matplotlib.pyplot as plt
import photutils
from photutils import DAOStarFinder
import pandas as pd
import pickle
import configparser
import glob

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
config.read("astrom/config.ini")

def find_stars(dateString,number_of_dithers):

    star_coords_every_dither = {} # initialize dictionary to hold the coords of found stars

    # initialize dictionary to contain x,y of stars in pixel space
    dictStars = {}

    # obtain the list of files which have been dewarped and de-rotated
    # according to the PA in the FITS headers
    asterism_frames_directory_retrieve = str(config["data_dirs"]["DIR_ASTERISM_DEROT"])
    asterism_frames_post_derot_names = list(glob.glob(os.path.join(asterism_frames_directory_retrieve, "*.fits")))

    for ditherPos in range(0,len(asterism_frames_post_derot_names)):
        print("------------------------------")
        print("Finding star positions in dither position "+str(ditherPos)+" ...")

        # read in image and header
        imageMedian, header = fits.getdata(asterism_frames_post_derot_names[ditherPos],
                                           0, header=True)
        
        # find star locations; input parameters may need some fine-tuning to get good results
        coordsAsterism = find_pinhole_centroids.find_psf_centers(imageMedian,20.,800.)
        coordsAsterismBright = find_pinhole_centroids.find_psf_centers(imageMedian,20.,3000.) # this is just for guiding the eye in identifying the stars
        
        if (len(coordsAsterism)>1): # if there is >1 star, there are >=1 baselines

            '''
            import ipdb; ipdb.set_trace()
            star_coords_every_dither[keyName] = np.transpose([coordsAsterism[:,0], coordsAsterism[:,1]])
            '''
            indicesArray = np.arange(np.size(coordsAsterism[:,0]))
            print("------------------------------")
            print("Please check the plot and note the indices from among the printed star positions which are TRUE "
                  "stars and are also IDENTIFIABLE (i.e., their RA and DEC is known. Then close the plot to continue.")
            print("------------------------------")
            print('Index, x, y:')
            dataFrame = pd.DataFrame(coordsAsterism, columns=['[x]','[y]'])
            print(dataFrame)
            #print(np.transpose([indicesArray,coordsAsterism[:,0],coordsAsterism[:,1]]))
        
            # plots for cross-checking found centroids with true stars (save plots manually as desired)
            plt.imshow(imageMedian, origin="lower", cmap="gray", vmin=0, vmax=1500)
            plt.scatter(coordsAsterism[:,0], coordsAsterism[:,1], color="yellow", label="Questionable detections")
            #import ipdb; ipdb.set_trace()
            # put index numbers next to star
            [plt.annotate(str(dataFrame.index.values[i]), (coordsAsterism[i,0], coordsAsterism[i,1]), fontsize=20, color="white") for i in range(len(dataFrame.index.values))] 
            plt.scatter(coordsAsterismBright[:,0], coordsAsterismBright[:,1], s=60, color="orange", label="Strong detections")
            plt.title("LMIRCam Trapezium observation, "+dateString+
                      "\n"+str(os.path.basename(asterism_frames_post_derot_names[ditherPos])))
            plt.legend()
            plt.show()

            # which coordinates should be pickled?
            print("------------------------------")
            print('Type the indices, separated by spaces, of the elements that are good')
            goodPts = [int(x) for x in input().split()] # user inputs
            print("------------------------------")
            print("Type the star names (in the above order, separated by spaces) as per the naming convention in star_coords.csv")
            namesGoodPts = [str(x) for x in input().split()] # user inputs
            print("------------------------------")

            # append the coordinate info
            keyName = "dither_pos_"+"%02i"%ditherPos # name of dictionary key for this dither
            dictStars[keyName] = [goodPts, namesGoodPts, dataFrame]

        else:
            print('Did not find >=1 stellar pairs!') # number of baselines is zero
            print("------------------------------")
            
    # pickle
    pickle.dump(dictStars, open("identified_stars_"+dateString+".p", "wb"))
    print("Star dictionary pickled as "+"identified_stars_"+dateString+".p")
    print("------------------------------")
