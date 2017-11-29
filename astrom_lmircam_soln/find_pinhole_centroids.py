# This find the centers of the pinhole PSFs by putting down a grid of initial guesses and finding
# nearby maxima

# created 22 Nov 2016 by E.S.

import numpy as np
import scipy
from astropy.io import fits
import matplotlib.pyplot as plt
import photutils
from photutils import DAOStarFinder

def rotation_matrix(angleDegPass): # use degrees in argument
    radEquiv = angleDegPass*np.pi/180.
    rot = np.matrix([[np.cos(radEquiv), -np.sin(radEquiv)], \
                    [np.sin(radEquiv), np.cos(radEquiv)]])

    return rot

def put_down_grid_guesses(approxHoleSpacingPass,rotationAnglePass):
    # puts down a simple grid that puts points as closely as possible to the pinhole images
    xHoles = np.arange(9.0, 1499.0, approxHoleSpacingPass)
    yHoles = np.arange(165.0, 1655.0, approxHoleSpacingPass)
    xHolesMeshGrid, yHolesMeshGrid = np.meshgrid(xHoles,yHoles)
    coordMatrix = np.matrix(np.transpose([np.ravel(xHolesMeshGrid),np.ravel(yHolesMeshGrid)])) # rotate

    xx = (np.matmul(coordMatrix,rotation_matrix(rotationAnglePass)))[:,0]
    yy = (np.matmul(coordMatrix,rotation_matrix(rotationAnglePass)))[:,1]
    
    return xx, yy

def find_psf_centers(imagePass,fwhmPass,thresholdPass):
    # uses the grid of guesses to find the real centroids
    
    daofind = DAOStarFinder(fwhm=fwhmPass, threshold=thresholdPass, exclude_border=True)
    sources = daofind(imagePass)

    return sources['xcentroid'], sources['ycentroid']

def consistency_xy_found_guesses(xCoordsGuessesPass,yCoordsGuessesPass,xCoordsFoundPass,yCoordsFoundPass):
    # makes each of the entries in the x,y lists of the guessed and found pinholes refer to the same pinholes
    xCoordsGuesses_finalSortedPass = [] # initialize
    yCoordsGuesses_finalSortedPass = []
    xCoordsFound_finalSortedPass = []
    yCoordsFound_finalSortedPass = []

    for p in range(0,32): # loop over stripes in the pinhole image
        colSpacingPixels = 48
        indicesStripeGuesses = np.where((xCoordsGuessesPass > p*colSpacingPixels) & (xCoordsGuessesPass < (p+1)*colSpacingPixels))
        indicesStripeFound = np.where((xCoordsFoundPass > p*colSpacingPixels) & (xCoordsFoundPass < (p+1)*colSpacingPixels))
        # make subset of points that are within the stripe
        xCoordsGuessesSubset = xCoordsGuessesPass[indicesStripeGuesses]
        yCoordsGuessesSubset = yCoordsGuessesPass[indicesStripeGuesses]
        xCoordsFoundSubset = xCoordsFoundPass[indicesStripeFound]
        yCoordsFoundSubset = yCoordsFoundPass[indicesStripeFound]
        # sort subsets of points in y, and tack onto arrays
        xCoordsGuessesSubset_sortedinY = (np.ravel(xCoordsGuessesSubset))[np.ravel(np.argsort(yCoordsGuessesSubset))]
        yCoordsGuessesSubset_sortedinY = (np.ravel(yCoordsGuessesSubset))[np.ravel(np.argsort(yCoordsGuessesSubset))]
        xCoordsFoundSubset_sortedinY = (np.ravel(xCoordsFoundSubset))[np.ravel(np.argsort(yCoordsFoundSubset))]
        yCoordsFoundSubset_sortedinY = (np.ravel(yCoordsFoundSubset))[np.ravel(np.argsort(yCoordsFoundSubset))]
        # concatenate
        xCoordsGuesses_finalSortedPass = np.concatenate((xCoordsGuesses_finalSortedPass,xCoordsGuessesSubset_sortedinY))
        yCoordsGuesses_finalSortedPass = np.concatenate((yCoordsGuesses_finalSortedPass,yCoordsGuessesSubset_sortedinY))
        xCoordsFound_finalSortedPass = np.concatenate((xCoordsFound_finalSortedPass,xCoordsFoundSubset_sortedinY))
        yCoordsFound_finalSortedPass = np.concatenate((yCoordsFound_finalSortedPass,yCoordsFoundSubset_sortedinY))

    return xCoordsGuesses_finalSortedPass, yCoordsGuesses_finalSortedPass,\
        xCoordsFound_finalSortedPass,yCoordsFound_finalSortedPass
