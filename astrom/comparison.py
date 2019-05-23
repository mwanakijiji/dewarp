# This finds star locations in pixel space in asterism images that have been dewarped (using the
# solution found by find_dewarp_solution.py) and derotated (based on the PA in the FITS header).
# The coordinates (including false positives) are printed to the Terminal and have to be manually
# cross-checked with known star positions, after which the pixel coordinates have to be written
# into the appropriate dictionaries in find_plate_scale_and_orientation.py.

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import pandas as pd
import datetime
from astropy.io import fits
import asciitable
from pathlib import Path
import pickle
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
import math
import itertools

def cdf_fcn(array_input):
    '''
    Return CDF of an unsorted input array of values
    '''
    
    number_cum_norm = np.divide(np.arange(len(array_input)),len(array_input))
    array_input_sort = np.sort(array_input)
    array_cdf = np.divide(np.cumsum(array_input_sort),np.cumsum(array_input_sort)[-1])

    return array_input_sort, number_cum_norm


def order_Y(star1input, star2input, star1nameInput, star2nameInput):
    '''
    determine which object is higher in y in (x,y) space
    '''
    if (star1input[1] > star2input[1]):
        c_high = star1input
        c_low = star2input
        name_high = star1nameInput
        name_low = star2nameInput
    else:
        c_high = star2input
        c_low = star1input
        name_high = star2nameInput
        name_low = star1nameInput
        
    # find position angle (in deg E of N)
    del_x = -np.subtract(c_high[0],c_low[0]) # minus sign because we want del_x E of N
    del_y = np.subtract(c_high[1],c_low[1])
    pos_angle_XY = np.arctan(np.divide(del_x,del_y))*(180./np.pi)
        
    return c_low, c_high, pos_angle_XY, name_low, name_high


def angOffset_plateScale(dateString,plotTitle,plot=True):

    # read in pickled star positions in (x,y)
    picklefile = "identified_stars_"+dateString+".p"
    fo1=open(picklefile,'rb')
    dat=pickle.load(fo1)
    fo1.close()

    # read in text file of star positions in (RA,DEC)
    df = pd.read_csv('star_ra_dec_list.csv')

    # initialize quantities
    baselineNumber = 0 # for counting number of baselines
    angleRadecMinusXYArray = [] # for collecting angle differences
    plateScaleArray = [] # for collecting plate scale measurements

    # loop over dither positions
    for ditherPos in range(0,len(dat.keys())):
        
        keyName = "dither_pos_"+"%02i"%ditherPos # key for this dither
        parent = dat[keyName][1] # N star names

        # list of all combinations of names (no degeneracies)
        allCombs = list(itertools.combinations(parent, 2))
    
        # loop over all baselines (i.e., take N, pick 2)
        for baseline in range(0,len(allCombs)):
        
        
            ## find position angle and distance in (x,y) space
        
            # retrieve x, y
            star1name = allCombs[baseline][0]
            star2name = allCombs[baseline][1]
        
            names = np.array(dat[keyName][1][:])
            star1elem = int(np.where(names==star1name)[0]) # element number in names array
            star2elem = int(np.where(names==star2name)[0])
        
            star1elem_2 = dat[keyName][0][star1elem]
            star2elem_2 = dat[keyName][0][star2elem]
        
            star1coords_xy = [dat[keyName][2]['[x]'][star1elem_2],dat[keyName][2]['[y]'][star1elem_2]]
            star2coords_xy = [dat[keyName][2]['[x]'][star2elem_2],dat[keyName][2]['[y]'][star2elem_2]]

            # sort according to y-position
            coords_low, coords_high, pos_angle_xy, name_low, name_high = order_Y(
                star1coords_xy, star2coords_xy, star1name, star2name)
        
            # find distance between the two stars
            del_y = np.subtract(coords_high[1],coords_low[1])
            del_x = np.subtract(coords_high[0],coords_low[0])
            dist_xy = np.sqrt(np.power(del_x,2)+np.power(del_y,2))
        
        
            ## find position angle in (RA,DEC) space
            # (note that stars are in same order, which is important if they have nearly equal DEC or y)
        
            # retrieve RA, DEC
            radecStarElem_low = np.where(df[' shorthand']==' '+name_low)[0] # find element number
            radecStarElem_high = np.where(df[' shorthand']==' '+name_high)[0]
            raString_low = df[' RA'][radecStarElem_low].values
            decString_low = df[' DEC'][radecStarElem_low].values
            raString_high = df[' RA'][radecStarElem_high].values
            decString_high = df[' DEC'][radecStarElem_high].values

            # find angle, separation
            c_low = SkyCoord(raString_low+decString_low, unit=(u.hourangle, u.deg))
            c_high = SkyCoord(raString_high+decString_high, unit=(u.hourangle, u.deg))
            pos_angle_radec = c_low.position_angle(c_high).degree[0] # position angle, E of N
            sep_radec = c_low.separation(c_high).arcsec[0] # separation in asec
        
            baselineNumber += 1 # chalk up this baseline to the total
        
            # how much further E of W is the position angle from (RA,DEC) than (x,y)?
            angleDiff_1 = np.subtract(pos_angle_radec,pos_angle_xy)
        
            if (angleDiff_1 > 0): # if x,y angle opens further E of N than the RA,DEC angle
                angleDiff = np.mod(
                    angleDiff_1,
                    360.) # mod is in case one angle is <0 and the other >180
            else: # if difference between x,y angle and RA,DEC angle is negative
                angleDiff = np.copy(angleDiff_1)
            
            angleRadecMinusXYArray = np.append(angleRadecMinusXYArray,angleDiff) # append del_angle to array
            plateScaleArray = np.append(plateScaleArray,1000.*np.divide(sep_radec,dist_xy)) # append plate scale (mas/pix)

            
    ## make CDFs...
    # ...of angular offsets
    #angleDiffArraySorted = sorted(angleRadecMinusXYArray)
    #angleDiff_csf = np.cumsum(angleDiffArraySorted).astype("float32")
    #angleDiff_csf_norm = np.divide(angleDiff_csf,np.max(angleDiff_csf))

    angleDiffArraySorted, number_angleDiff_cum_norm = cdf_fcn(angleRadecMinusXYArray)
    
    # ...of plate scales
    #plateScaleArraySorted = sorted(plateScaleArray)
    #plateScale_csf = np.cumsum(plateScaleArraySorted).astype("float32")
    #plateScale_csf_norm = np.divide(plateScale_csf,np.max(plateScale_csf))

    plateScaleArraySorted, number_plateScale_cum_norm = cdf_fcn(plateScaleArray)

    ## find median, +- sigma values...
    # ...of angular offsets
    angleDiff_negSigmaPercentile = np.percentile(angleRadecMinusXYArray,15.9)
    angleDiff_50Percentile = np.percentile(angleRadecMinusXYArray,50)
    angleDiff_posSigmaPercentile = np.percentile(angleRadecMinusXYArray,84.1)
    # ...of plate scales
    plateScale_negSigmaPercentile = np.percentile(plateScaleArray,15.9)
    plateScale_50Percentile = np.percentile(plateScaleArray,50)
    plateScale_posSigmaPercentile = np.percentile(plateScaleArray,84.1)

    # prepare strings
    string1 = '{0:.3f}'.format(angleDiff_50Percentile)
    string2 = '{0:.3f}'.format(np.subtract(angleDiff_posSigmaPercentile,angleDiff_50Percentile))
    string3 = '{0:.3f}'.format(np.subtract(angleDiff_50Percentile,angleDiff_negSigmaPercentile))
    string4 = '{0:.3f}'.format(plateScale_50Percentile)
    string5 = '{0:.3f}'.format(np.subtract(plateScale_posSigmaPercentile,plateScale_50Percentile))
    string6 = '{0:.3f}'.format(np.subtract(plateScale_50Percentile,plateScale_negSigmaPercentile))

    # print info
    print('------------------------------')
    print('Number of stellar pair baselines:')
    print(baselineNumber)
    print('------------------------------')
    print('Need to rotate array E of N:\n'+string1+'/+'+string2+'/-'+string3+' deg')
    print('------------------------------')
    print('Plate scale:')
    print(string4+'/+'+string5+'/-'+string6+' mas/pix')

        
    if (plot):
    
        # plot rotation angle
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.axvline(x=angleDiff_negSigmaPercentile,linestyle='--',color='k')
        ax.axvline(x=angleDiff_50Percentile,linestyle='-',color='k')
        ax.axvline(x=angleDiff_posSigmaPercentile,linestyle='--',color='k')
        ax.scatter(angleDiffArraySorted, number_angleDiff_cum_norm)
        ax.text(0.8, 0.05,s='Need to rotate array E of N:\n'+string1+'/+'+string2+'/-'+string3+' deg\n\nStellar baselines:\n'+str(baselineNumber))
        plt.title('CDF of difference (E of N) between (RA, DEC) and (x, y) position angles on LMIRcam, '+plotTitle)
        plt.xlabel('Degrees E of N')
        plt.ylabel('Normalized CDF')
        plt.show()

        # plot plate scale
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.axvline(x=plateScale_negSigmaPercentile,linestyle='--',color='k')
        ax.axvline(x=plateScale_50Percentile,linestyle='-',color='k')
        ax.axvline(x=plateScale_posSigmaPercentile,linestyle='--',color='k')
        ax.scatter(plateScaleArraySorted, number_plateScale_cum_norm)
        ax.text(10.4, 0.8,s='Plate scale:\n'+string4+'/+'+string5+'/-'+string6+' mas/pix\n\nStellar baselines:\n'+str(baselineNumber))
        plt.title('Plate scale of LMIRcam, '+plotTitle)
        plt.xlabel('PS (mas/pix)')
        plt.ylabel('Normalized CDF')
        plt.show()
