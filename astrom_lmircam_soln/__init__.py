# created by E.S., 22 Nov 2016
# revised 2017 Dec

import os
import subprocess
import numpy as np
from scipy.ndimage import rotate, map_coordinates

path=os.path.dirname(__file__)

dateStringAsterism = 'ut_2017_11_06'
dateStringPinholes = 'ut_2017_11_08'
generalStem = os.path.expanduser('~')+'/../../media/unasemaje/Seagate Expansion Drive/lbti_data_reduction/lmircam_astrometry/'

# data paths: pinholes
'''
raw_pinholes_data_stem = generalStem+dateStringPinholes+'/pinholeGrid/rawData/' # obtain raw science data only (don't save anything to this!)
calibrated_pinholes_data_stem = "dx_pinholes_190125.fits" # deposit science arrays after bias-subtraction, flat-fielding, etc.
#save_pinholes_data_stem = generalStem+'/textfile_results/' # deposit polynomial fit data

# data paths: trapezium
raw_trapezium_data_stem = generalStem+dateStringAsterism+'/asterism/rawData/' # obtain raw science data only (don't save anything to this!)
calibrated_trapezium_data_stem = generalStem+dateStringAsterism+'/asterism/processedData/' # deposit science arrays after bias-subtraction, flat-fielding, etc.
#save_trapezium_data_stem = generalStem+'/textfile_results/' # deposit polynomial fit data
'''


def get_git_hash():
    '''
    returns the git hash for recording the exact version of the pipeline used
    '''
    gitd=path+'/../' # directory containing .git
    hash=subprocess.check_output('git --git-dir=%s/.git --work-tree=%s rev-parse HEAD'%(gitd,gitd),shell=True)
    return hash.strip()


def write_hash(filenameString):
    '''
    write the git hash
    '''
    fo=open(filenameString,'w')
    fo.write(str(get_git_hash()))
    fo.write('\n')
    fo.close()


# this function from Jordan Stone's LEECH code
def rot(im, angle, axis, order=3, pivot=False):
    '''rotate an image clockwise by angle [degrees] about axis.
    if pivot is true the image will pivot about the axis. otherwise
    the axis will be centered in the output image''' 
    angle*=np.pi/180.#convert to radians
    y,x = np.indices(im.shape)

    # calculate how the axis moves when pivoting from bottom left corner
    theta_axis = np.arctan2(axis[1],axis[0])
    r_axis = np.abs(axis[0]+1j*axis[1])
    yoffset = r_axis*np.sin(theta_axis) - r_axis*np.sin(theta_axis-angle)
    xoffset = r_axis*np.cos(theta_axis) - r_axis*np.cos(theta_axis-angle)

    # put the axis in the middle? 
    ycenter_offset = (1-pivot) * ((im.shape[0]/2.)-axis[1])#pivot is a bool (i.e. 0 or 1)
    xcenter_offset = (1-pivot) * ((im.shape[1]/2.)-axis[0])
    yoffset += ycenter_offset
    xoffset += xcenter_offset

    # make rotation matrix elements
    ct = np.cos(angle)
    st = np.sin(angle)

    # do the rotation 
    new_x = (ct*(x-xoffset) - st*(y-yoffset)) 
    new_y = (st*(x-xoffset) + ct*(y-yoffset))
    
    return map_coordinates(im, [new_y, new_x], order=order)
