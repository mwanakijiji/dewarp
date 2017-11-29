import numpy as np
from scipy.ndimage import map_coordinates
import scipy.spatial as spatial
import scipy.spatial.distance as dist

# set down grid to sample the detector area
def barb_grid():
    # puts down a grid to sample the detector area; very similar to put_down_grid_guesses, except that spacing and rotation (none!) are hardwired
    spacing = 48.0 # this many pixels between the holes
    #rotationAnglePass = 0.0
    xHoles = np.arange(0.5*spacing, 43.5*spacing, spacing)
    yHoles = np.arange(0.5*spacing, 43.5*spacing, spacing) #43.5*spacing, spacing) ## CHANGE MADE
    xHolesMeshGrid, yHolesMeshGrid = np.meshgrid(xHoles,yHoles)
    #coordMatrix = np.matrix(np.transpose([np.ravel(xHolesMeshGrid),np.ravel(yHolesMeshGrid)]))

    #xx = coordMatrix[:,0]
    #yy = coordMatrix[:,1]
    
    return xHolesMeshGrid, yHolesMeshGrid, xHoles, yHoles


# separate barb locations into those that are close to, and far away from, the region sampled by the pinhole grid
def find_close_pts_cdist(points1, points2, cutoff):
    
    indicesAll = np.where(dist.cdist(points1, points2))[0]
    indicesClose = np.where(dist.cdist(points1, points2) <= cutoff)[0]
    indicesClose = np.unique(indicesClose)
    indicesFar = indicesAll[np.in1d(indicesAll,indicesClose,invert=True)]
    indicesFar = np.unique(indicesFar)
    
    return indicesClose, indicesFar
