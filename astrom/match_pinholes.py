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
import configparser

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
config.read("astrom/config.ini")



#####################################################################
# RETRIEVE MEDIAN PINHOLE GRID IMAGE, FIND PINHOLE CENTERS

def match_pinholes(translationPass,
                   holeSpacingPass,
                   barrelCenterPass,
                   barrelAmountPass,
                   rotationAnglePass,
                   writeoutString='test',
                   plotTitleString='test',
                   plot=True):
    '''
    translationPass: translation of the grid
    holeSpacingPass: spacing between the holes
    barrelCenterPass: TBD
    barrelCenterPass: TBD
    barrelCenterPass: TBD
    writeoutString: TBD
    writeoutString: TBD
    '''            

    # write the git hash
    write_hash('git_hash_match_pinholes.txt')

    # read in median image of pinholes
    hdul = fits.open(config["data_dirs"]["DIR_PINHOLE_BASIC"] + config["src_file_names"]["PINHOLE_FITS"]) 
    imagePinholes = hdul[0].data.copy()
    
    # make model pinhole locations: distorted and undistorted
    coordsModel_d, coordsModel_not_d = find_pinhole_centroids.put_down_grid_guesses(
        translationPass,
        holeSpacingPass,
        barrelCenterPass,
        barrelAmountPass,
        rotationAnglePass)

    # find empirical pinhole locations
    coordsEmpirical = find_pinhole_centroids.find_psf_centers(
        imagePass = imagePinholes,
        fwhmPass = 20.,
        thresholdPass = 1000.)

    # option for adding in some additional points manually
    ##xCoordsMissed = [71.774,697.353,1460.66]
    ##yCoordsMissed = [1267.57,1404.06,737.932]
    '''
    xCoordsFound = np.concatenate((xCoordsFoundAutomated,
                               xCoordsMissed),
                              axis=0)
    yCoordsFound = np.concatenate((yCoordsFoundAutomated,
                               yCoordsMissed),
                              axis=0)
    '''

    # match and sort the model coordinates (distorted and undistorted) with the empirical pinhole coordinates
    coordsModel_d_matched, coordsModel_not_d_matched, coordsEmpirical_matched = find_pinhole_centroids.match_model_empirical(
        coordsModel_d,
        coordsModel_not_d,
        coordsEmpirical,
        imagePinholes,
        plotTitleString,
        plot=plot)

    # pickle
    picklename='pinhole_coordinates_'+writeoutString+'.pkl'
    fo=open(picklename,'wb')
    pickle.dump({'coordsModel_d_matched':coordsModel_d_matched,
                 'coordsModel_not_d_matched':coordsModel_not_d_matched, 
                 'coordsEmpirical_matched':coordsEmpirical_matched},
                 fo)
    fo.close()

    print('------------------------------')
    print('Saved pinhole grid info to:')
    print(picklename)
    print('------------------------------')

'''
class match:
    def __init__(self,translation, holeSpacing, barrelCenter, barrelAmount, rotationAngle, writeoutString, plot=True):
        self.translation = translation
        self.holeSpacing = holeSpacing
        self.barrelCenter = barrelCenter
        self.barrelAmoung = barrelAmount
        self.rotationAngle = rotationAngle
        self.writeoutString = writeoutString

    def __call__(self):
        match_pinholes(self.translation,
                       self.holeSpacing,
                       self.barrelCenter,
                       self.barrelAmount,
                       self.rotationAngle,
                       self.writeoutString)

match
'''
