# This find the centers of the pinhole PSFs by putting down a grid of initial guesses and finding
# nearby maxima

# created 22 Nov 2016 by E.S.

import numpy as np
import scipy
from astropy.io import fits
import matplotlib.pyplot as plt
import photutils
from photutils import DAOStarFinder
from scipy import spatial

def rotation_matrix(angleDegPass): # use degrees in argument
    radEquiv = angleDegPass*np.pi/180.
    rot = np.matrix([[np.cos(radEquiv), -np.sin(radEquiv)], \
                    [np.sin(radEquiv), np.cos(radEquiv)]])

    return rot

def put_down_grid_guesses(displacementPass,approxHoleSpacingPass,barrelCenterPass,barrelDegreePass,rotationAnglePass):
    '''
    puts down a simple grid that puts points as closely as possible to the pinhole images
    INPUTS:
    displacementPass: amount to translate the model grid ([x,y], in pixels)
    approxHoleSpacingPass: the inter-pinhole spacing (in pixels)
    barrelCenterPass: the center of the barrel distortion in the unrotated frame ([x,y], in pixels)
    barrelDegreePass: the amount of barrel distortion
    rotationAnglePass: how much does the ideal grid have to be rotated to overlay it on the observed pinholes? (in deg)
    '''
    #xHoles = np.arange(9.0, 1499.0, approxHoleSpacingPass)
    #yHoles = np.arange(165.0, 1655.0, approxHoleSpacingPass)
    xHoles = np.arange(0.0, 2048.0, approxHoleSpacingPass)
    yHoles = np.arange(0.0, 2048.0, approxHoleSpacingPass)
    xHolesMeshGrid, yHolesMeshGrid = np.meshgrid(xHoles,yHoles) # arrange into an (x,y) grid
    xHolesMeshGrid_no_d = np.copy(xHolesMeshGrid) # make copy for undistorted coordinates
    yHolesMeshGrid_no_d = np.copy(yHolesMeshGrid)
    
    # add barrel distortion: radial distance from a point (x0,y0) is r_d = r_u*(1-alpha*|r_u|**2) 
    xHolesMeshGrid_d, yHolesMeshGrid_d = np.multiply(
        [xHolesMeshGrid-barrelCenterPass[0], yHolesMeshGrid-barrelCenterPass[1]], # x_u-x_0, y_u-y_0
        np.subtract(
            1.,
            np.multiply(
                np.power( np.subtract(xHolesMeshGrid,barrelCenterPass[0]), 2 )  +  np.power( np.subtract(yHolesMeshGrid,barrelCenterPass[1]), 2 ),
                barrelDegreePass
            )
        )
        )

    # add (x0,y0) back in, as well as the desired translation
    xHolesMeshGrid_d += barrelCenterPass[0] + displacementPass[0]
    yHolesMeshGrid_d += barrelCenterPass[1] + displacementPass[1]
    xHolesMeshGrid_no_d += displacementPass[0] # note barrelCenterPass was not subtracted out
    yHolesMeshGrid_no_d += displacementPass[1]

    # rotate
    coordMatrix_d = np.matrix(np.transpose([np.ravel(xHolesMeshGrid_d),np.ravel(yHolesMeshGrid_d)]))
    coordMatrix_no_d = np.matrix(np.transpose([np.ravel(xHolesMeshGrid_no_d),np.ravel(yHolesMeshGrid_no_d)]))
    xHolesMeshGrid_d_r = (np.matmul(coordMatrix_d,rotation_matrix(rotationAnglePass)))[:,0]
    yHolesMeshGrid_d_r = (np.matmul(coordMatrix_d,rotation_matrix(rotationAnglePass)))[:,1]
    xHolesMeshGrid_no_d_r = (np.matmul(coordMatrix_no_d,rotation_matrix(rotationAnglePass)))[:,0]
    yHolesMeshGrid_no_d_r = (np.matmul(coordMatrix_no_d,rotation_matrix(rotationAnglePass)))[:,1]

    ptsModel_d_Pass = np.concatenate((xHolesMeshGrid_d_r,yHolesMeshGrid_d_r),axis=1)
    ptsModel_no_d_Pass = np.concatenate((xHolesMeshGrid_no_d_r,yHolesMeshGrid_no_d_r),axis=1)

    # return distorted [x,y], undistorted [x,y]
    return ptsModel_d_Pass, ptsModel_no_d_Pass

def find_psf_centers(imagePass,fwhmPass,thresholdPass):
    '''
    employ DAOPhot to find PSF coordinates
    '''
    daofind = DAOStarFinder(fwhm=fwhmPass, threshold=thresholdPass, exclude_border=True)
    sources = daofind(imagePass)

    # put the x,y coordinates into an array
    ptsEmpiricalPass = np.transpose(np.concatenate(([sources['xcentroid']], [sources['ycentroid']]),axis=0))

    return ptsEmpiricalPass


def match_model_empirical(
        ptsModel_d_Pass,
        ptsModel_not_d_Pass,
        ptsEmpiricalPass,
        imagePinholesPass,
        plot=True):
    '''
    matches nearest neighbors between model and empirical pinhole locations
    INPUTS:
    ptsModel_d_Pass: the pinhole locations in the distorted model; needed for matching empirical pinholes on a nearest-neighbor basis (doesn't need to be sorted, but should be larger than ptsEmpiricalPass to make use of as many pinholes as possible)
    ptsModel_not_d_Pass: same as ptsModel_d_Pass, except free of distortion; needed to underpin dewarp solution
    ptsEmpiricalPass: the pinhole locations found by DAOPhot (doesn't need to be sorted)
    plot: display a plot as a check?
    '''

    treeModel = spatial.KDTree(ptsModel_d_Pass) # set up tree to query
    indicesOfInterest = treeModel.query(ptsEmpiricalPass) # find nearest neighbors; indicesOfInterest[1][:] are array indices
    xEmpirical = ptsEmpiricalPass[:,0] # rename
    yEmpirical = ptsEmpiricalPass[:,1]
    xModel_d = ptsModel_d_Pass[indicesOfInterest[1][:]][:,0] # sort 
    yModel_d = ptsModel_d_Pass[indicesOfInterest[1][:]][:,1]
    xModel_not_d = ptsModel_not_d_Pass[indicesOfInterest[1][:]][:,0] # apply same sorted indices as before to undistorted model
    yModel_not_d = ptsModel_not_d_Pass[indicesOfInterest[1][:]][:,1]

    if plot: # make a plot to check mapping between coordinates
        
        # image background
        plt.imshow(imagePinholesPass, origin="lower", cmap="gray")
        
        # plot coordinates
        plt.scatter(xEmpirical,yEmpirical,color="red")
        plt.scatter(xModel_d,yModel_d,color="yellow") # distorted
        plt.scatter(xModel_not_d,yModel_not_d,color="blue") # not distorted
        
        # list comprehension for drawing white lines between model and empirical points
        [plt.plot([xEmpirical[j],xModel_not_d[j]],[yEmpirical[j],yModel_not_d[j]], color="w") for j in range(len(xEmpirical))]

        plt.title('red = empirical; yellow = distorted model; blue = undistorted model')
        plt.xlabel('x (pix)')
        plt.ylabel('y (pix)')
        plt.show()

    # package everything
    coords_d_Pass = np.transpose(np.squeeze([xModel_d,yModel_d]))
    coords_not_d_Pass = np.transpose(np.squeeze([xModel_not_d,yModel_not_d]))
    coordsEmpiricalPass = np.transpose([xEmpirical,yEmpirical])

    return coords_d_Pass, coords_not_d_Pass, coordsEmpiricalPass
    

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
