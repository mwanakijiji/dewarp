# This finds the position angles between different objects in LMIRCam data and
# compares them to 'true' angles in order to determine the directional orientation
# of the detector. The dictionaries in the code below have to be populated manually,
# using data in pixel space from find_asterism_star_locations.py, and in RA-DEC space
# from astrometric sources (like Close+ 2012 ApJ 749:180)

# created 14 Dec 2016 by E.S.

import numpy as np
import scipy
from astropy.io import fits
import matplotlib.pyplot as plt
import photutils
from photutils import DAOStarFinder

# note all the found centroids (true and false positives) have been pickled in a binary file,
# in case you want to check it out manually


######################################################
# DICTIONARIES WITH STAR POSITIONS IN PIXEL SPACE

# note that these coordinates are from images that have been
# dewarped and derotated

detector_star_pos_dither_pos_0 = {
    "beta": [ 1278.77543474, 1320.9970804 ],
    "gamma": [ 1430.93666334, 1322.57580845 ], 
    "delta": [ 1273.55056779, 1029.67652732 ],
    "epsilon": [ 1349.51822334, 938.26286875 ],
    "sigma": [ 34.32329057, 1035.99366435 ]
}


detector_star_pos_dither_pos_1 = {
    "beta": [ 1494.83417109, 1158.32330349 ],
    "delta": [ 1489.15063578, 867.05557534 ],
    "epsilon": [ 1564.57084883, 775.90279951 ],
    "sigma": [ 250.39790634, 873.1764211 ]
}


detector_star_pos_dither_pos_2 = {} # no useful baselines in this dither


detector_star_pos_dither_pos_3 = {
    "zeta": [ 646.0806408, 1394.03423163 ],
    "eta": [ 831.16669266, 1424.42426591 ],
    "sigma": [ 690.31178058, 567.42001907 ],
    "theta": [ 1294.84738519, 1396.59477603 ]
}


detector_star_pos_dither_pos_4 = {
    "D1": [ 161.80662043, 1219.07114255 ],
    "zeta": [ 871.04551665, 1246.11397877 ],
    "eta": [ 1055.49264695, 1278.34187839 ],
    "theta": [ 1520.53672728, 1249.9000294 ],
    "sigma": [ 915.50243811, 420.80227143 ]
}


detector_star_pos_dither_pos_5 = {
    "D1": [ 388.06814022, 1078.17689859 ],
    "zeta": [ 1097.01870031, 1105.44197087 ],
    "eta": [ 1281.24396884, 1136.84787458 ],
    "theta": [ 1747.36109572, 1109.48549406 ],
    "sigma": [ 1142.77871476, 280.20717649 ]
}


detector_star_pos_dither_pos_6 = {
    "D1": [ 623.05894952, 945.58783417 ],
    "sigma": [ 1378.0252694, 149.03119895 ]
}


detector_star_pos_dither_pos_7 = {
    "D1": [ 260.82918971, 863.17440254 ],
    "zeta": [ 970.22724374, 891.19689695 ],
    "theta": [ 1619.54197816, 894.55645244 ],
    "sigma": [ 1017.48077575, 67.74322368 ]
}


detector_star_pos_dither_pos_8 = {
    "D1": [ 268.39410695, 866.42409421 ],
    "theta": [ 1626.62154931, 898.05974716 ]
}


detector_star_pos_dither_pos_9 = {
    "B1": [ 1468.3469216, 1729.63239681 ],
    "B2": [ 1554.27095054, 1707.74522771 ],
    "zeta": [ 619.28447452, 802.77405726 ],
    "eta": [ 802.37654042, 833.40172288 ],
    "theta": [ 1267.61175361, 805.3130301 ]
}


detector_star_pos_dither_pos_10 = {
    "B1": [ 1109.60884172, 1631.33099305 ],
    "B2": [ 1195.36298066, 1609.53788767 ],
    "zeta": [ 263.09348992, 704.86976588 ],
    "eta": [ 447.43134426, 736.14584106 ],
    "theta": [ 911.37798525, 707.22974711 ],
    "E1": [ 1622.75226422, 1349.94401585 ],
    "A1": [ 1556.35404063, 929.41360769]
}


detector_star_pos_dither_pos_11 = {
    "B1": [ 752.41186969, 1526.55925515 ],
    "B2": [ 838.11725978, 1504.50337841 ],
    "beta": [ 1194.64599105, 61.37835351 ],
    "theta": [ 554.36176368, 603.35113074 ],
    "E1": [ 1265.22196811, 1245.174892 ],
    "A1": [ 1198.64372109, 824.78809778 ]
}


detector_star_pos_dither_pos_12 = {
    "B1": [ 627.9749373, 1287.25096115 ],
    "B2": [ 713.98189207, 1264.87224362 ],
    "theta": [ 431.58300764, 364.45026387 ],
    "E1": [ 1140.66986209, 1005.94140448 ],
    "A1": [ 1075.2046255, 585.79380025 ]
}


detector_star_pos_dither_pos_13 = {
    "B1": [ 505.51340081, 1042.9600791 ],
    "B2": [ 591.27823717, 1021.12844839 ],
    "theta": [ 308.31797276, 120.49831303 ],
    "E1": [ 1018.15728133, 761.72791744 ],
    "A1": [ 952.59468061, 341.54887559 ]
}


detector_star_pos_dither_pos_14 = {
    "B1": [ 388.75557516, 802.82127489 ],
    "B2": [ 474.41997228, 780.90329529 ],
    "E1": [ 901.41867217, 522.01639754 ],
    "A1": [ 835.85676622, 101.6476346 ]
}


detector_star_pos_dither_pos_15 = {
    "B1": [ 272.10697931, 560.18210666 ],
    "B2": [ 358.38133726, 538.27847944 ],
    "E1": [ 785.45682778, 279.2418408 ]
}


detector_star_pos_dither_pos_16 = {
    "B2": [ 257.59547847, 291.42699017 ],
    "E1": [ 685.11940769, 32.18144383 ]
}


detector_star_pos_dither_pos_17 = {} # no useful baselines in this dither


######################################################
# DICTIONARIES WITH STAR POSITIONS IN RA-DEC SPACE

# position vectors according to Close+ 2012 ApJ 749:180

# entries contain
# [0]: RA, in seconds (prefix '5:35:' to get full RA has been removed)
# [1]: DEC, in seconds of arc (prefix '-5:23:' to get full DEC has been removed)

true_star_pos = {
    "beta": [15.8337, 22.4207],
    "gamma": [15.7255, 22.4347],
    "delta": [15.8408, 25.5078],
    "epsilon": [15.7879, 26.5168],
    "zeta": [16.7469, 16.3777],
    "eta": [16.6148, 16.0836],
    "theta": [16.283, 16.512],
    "sigma": [16.7236, 25.1688],
    "B2": [16.069, 6.96452],
    "A1": [15.8202, 14.2891],
    "B1": [16.1299, 6.71895],
    "D2": [17.1675, 17.0013],
    "D1": [17.2558, 16.5298],
    "E1": [15.7673, 9.82764]
}


######################################################
# DEFINE FUNCTIONS

# -------------------------------------------------------------------- #
# fcn returns distance in RA, in sky-projected arcsecs
# INPUT:
# obj1Pass, obj2Pass: two dictionary entries for two objects
# flag 'lmir' should be set if the inputs are referring to objects in LMIRCam pixel space
# OUTPUT:
# distance in RA, in sky-projected arcsecs (the two objects can be entered in any
# order, but the output sign convention is such that the difference is
# del_RA = (RA_f - RA_i)*cos(DEC), where RA_f refers to the object furthest north)
def delta_RA(obj1Pass, obj2Pass, lmir=False):

    if not lmir:
        if (obj1Pass[1] < obj2Pass[1]): # choose object further N as 'final', that further S as 'initial'; note that DEC is positive S here
            objHigh = obj1Pass
            objLow = obj2Pass
        else:
            objHigh = obj2Pass
            objLow = obj1Pass
            
        ra_f = objHigh[0]*15. # factor of 15 to convert from sec LST to arcsec of RA
        dec_f = objHigh[1]
        ra_i = objLow[0]*15.
        dec_i = objLow[1]
        
        avgDec = 0.5*(dec_i+dec_f)
        degreesDECSouth = 5.+23./60.+avgDec/3600. # note this is 'positive south'
        diff = np.multiply((ra_f-ra_i),np.cos(degreesDECSouth*np.pi/180.))
        
    if lmir: # if in pixel space, where higher y means further N
        if (obj1Pass[1] > obj2Pass[1]):
            objHigh = obj1Pass
            objLow = obj2Pass
        else:
            objHigh = obj2Pass
            objLow = obj1Pass
            
        ra_f = objHigh[0]
        dec_f = objHigh[1]
        ra_i = objLow[0]
        dec_i = objLow[1]
        
        avgDec = 16.5 # close enough; this is the value for star 'theta'
        degreesDECSouth = 5.+23./60.+avgDec/3600. # note this is 'positive south'
        diff = np.multiply((ra_f-ra_i),np.cos(degreesDECSouth*np.pi/180.))
        diff = -diff # add sign change (RA goes + as pixels go - in x)

    return diff


# -------------------------------------------------------------------- #
# fcn returns distance in DEC
# INPUT:
# obj1Pass, obj2Pass: two dictionary entries for two objects
# flag 'lmir' should be set if the inputs are referring to objects in LMIRCam pixel space
# a flip of 180 deg in the angular difference when LMIRCam and the 'true' positions disagree on which object is further N
# OUTPUT:
# distance in DEC, in sky-projected arcsecs (the two objects can be entered in any
# order, but the output sign convention is such that the difference is
# del_DEC = DEC_f - DEC_i, where DEC_f refers to the object furthest north)
def delta_DEC(obj1Pass, obj2Pass, lmir=False):
            
    if lmir: # if in pixel space, where +y means further N
        if (obj1Pass[1] > obj2Pass[1]):
            objHigh = obj1Pass
            objLow = obj2Pass
        else:
            objHigh = obj2Pass
            objLow = obj1Pass

    if (not lmir): # if coords are 'true'
        if (obj1Pass[1] < obj2Pass[1]): # choose object further N as 'final', that further S as 'initial'; note that DEC is positive S here
            objHigh = obj1Pass
            objLow = obj2Pass
        else:
            objHigh = obj2Pass
            objLow = obj1Pass

    dec_f = objHigh[1]
    dec_i = objLow[1]

    if not lmir:
        diff = -(dec_f-dec_i) # negative sign since the DEC values here are 'positive S'
    else:
        diff = (dec_f-dec_i) # in pixel space, where +y means N

    return diff


# -------------------------------------------------------------------- #
# fcn returns distance between the two objects and the angle theta
# INPUT:
# star1Pass, star2Pass: two dictionary entries for two objects
# flag 'lmir' should be set if the inputs are referring to objects in LMIRCam pixel space
# OUTPUT:
# distancePass: distance between the two objects, in sky-projected arcsec
# thetaDegreePass: angle offset of the line connecting the two objects, relative to the
# meridian (convention here is 'positive angle E of N' (take care that the input coords are from an image
# that has already been derotated so PA=0)
def dist_and_theta(star1Pass, star2Pass, lmir=False):

    thetaRad = np.arctan(np.divide(delta_RA(star1Pass, star2Pass, lmir=lmir),
                                    delta_DEC(star1Pass, star2Pass, lmir=lmir)))
    thetaDegPass = thetaRad*180./np.pi
    distancePass = np.sqrt(delta_DEC(star1Pass, star2Pass, lmir=lmir)**2 + 
                           delta_RA(star1Pass, star2Pass, lmir=lmir)**2)

    import ipdb; ipdb.set_trace()
    return distancePass, thetaDegPass # asec OR pixels, deg E of N (at 0 deg PA)


# -------------------------------------------------------------------- #
# fcn returns the LMIRCam plate scale and angle offset for one pair of objects
# INPUT:
# lmirKeyObj1Pass, lmirKeyObj2Pass: dictionary entries for two objects on LMIRCam readout (ex. 'detector_star_pos_dither_pos_0["beta"]')
# trueSkyObj1Pass, trueSkyObj2Pass: same as above, from 'true' astrometry
# OUTPUT:
# [plate scale, angle offset] (units asec/pix and degrees)
def plate_scale_and_angles(lmirKeyObj1Pass, lmirKeyObj2Pass,
                           trueSkyObj1Pass, trueSkyObj2Pass):
    
    plateScalePass = np.divide( 
        dist_and_theta(
            trueSkyObj1Pass,
            trueSkyObj2Pass
        )[0], # distance in arcsec
        dist_and_theta(
            lmirKeyObj1Pass,
            lmirKeyObj2Pass,
            lmir=True
        )[0] # distance in pixels from LMIRCam readout array
        ) # return asec/pixel
    
    anglePass = np.subtract(
        dist_and_theta(
            lmirKeyObj1Pass,
            lmirKeyObj2Pass,
            lmir=True
        )[1], # baseline angle E of N on LMIRCam at PA=0
       dist_and_theta(
            trueSkyObj1Pass,
            trueSkyObj2Pass
        )[1] # true angle
        ) # return deg in 'CCW' direction
    
    return [plateScalePass, anglePass]


######################################################
# DICTIONARIES WITH DIFFERENCES BETWEEN VECTORS
# CONTAINING PIXEL-SPACE AND TRUE DISTANCES AND ANGLES


dither_pos_0_baselines = {
    
    "dp_0_del_beta_gamma":
    plate_scale_and_angles(detector_star_pos_dither_pos_0["beta"],
                           detector_star_pos_dither_pos_0["gamma"],
                           true_star_pos["beta"], 
                           true_star_pos["gamma"]),
    
    "dp_0_del_beta_delta":
    plate_scale_and_angles(detector_star_pos_dither_pos_0["beta"], 
                           detector_star_pos_dither_pos_0["delta"],
                           true_star_pos["beta"],
                           true_star_pos["delta"]),
    
    "dp_0_del_beta_epsilon": 
    plate_scale_and_angles(detector_star_pos_dither_pos_0["beta"], 
                           detector_star_pos_dither_pos_0["epsilon"],
                           true_star_pos["beta"],
                           true_star_pos["epsilon"]),
    
    "dp_0_del_beta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_0["beta"], 
                           detector_star_pos_dither_pos_0["sigma"],
                           true_star_pos["beta"],
                           true_star_pos["sigma"]),
#--------------------------------------------------------    
    "dp_0_del_gamma_delta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_0["gamma"], 
                           detector_star_pos_dither_pos_0["delta"],
                           true_star_pos["gamma"],
                           true_star_pos["delta"]),
    
    "dp_0_del_gamma_epsilon": 
    plate_scale_and_angles(detector_star_pos_dither_pos_0["gamma"], 
                           detector_star_pos_dither_pos_0["epsilon"],
                           true_star_pos["gamma"],
                           true_star_pos["epsilon"]),
    
    "dp_0_del_gamma_sigma":  
    plate_scale_and_angles(detector_star_pos_dither_pos_0["gamma"], 
                           detector_star_pos_dither_pos_0["sigma"],
                           true_star_pos["gamma"],
                           true_star_pos["sigma"]),
    
#--------------------------------------------------------
    "dp_0_del_delta_epsilon": 
    plate_scale_and_angles(detector_star_pos_dither_pos_0["delta"], 
                           detector_star_pos_dither_pos_0["epsilon"],
                           true_star_pos["delta"],
                           true_star_pos["epsilon"]),
    
    "dp_0_del_delta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_0["delta"], 
                           detector_star_pos_dither_pos_0["sigma"],
                           true_star_pos["delta"],
                           true_star_pos["sigma"]),
    
#--------------------------------------------------------
    "dp_0_del_epsilon_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_0["epsilon"], 
                           detector_star_pos_dither_pos_0["sigma"],
                           true_star_pos["epsilon"],
                           true_star_pos["sigma"])
    }


dither_pos_1_baselines = {
    
    "dp_1_del_beta_delta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_1["beta"], 
                           detector_star_pos_dither_pos_1["delta"],
                           true_star_pos["beta"],
                           true_star_pos["delta"]),
    
    "dp_1_del_beta_epsilon": 
    plate_scale_and_angles(detector_star_pos_dither_pos_1["beta"], 
                           detector_star_pos_dither_pos_1["epsilon"],
                           true_star_pos["beta"],
                           true_star_pos["epsilon"]),
    
    "dp_1_del_beta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_1["beta"], 
                           detector_star_pos_dither_pos_1["sigma"],
                           true_star_pos["beta"],
                           true_star_pos["sigma"]),
    
#--------------------------------------------------------
    "dp_1_del_delta_epsilon": 
    plate_scale_and_angles(detector_star_pos_dither_pos_1["delta"], 
                           detector_star_pos_dither_pos_1["epsilon"],
                           true_star_pos["delta"],
                           true_star_pos["epsilon"]),
    
    "dp_1_del_delta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_1["delta"], 
                           detector_star_pos_dither_pos_1["sigma"],
                           true_star_pos["delta"],
                           true_star_pos["sigma"]),
    
#--------------------------------------------------------
    "dp_1_del_epsilon_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_1["epsilon"], 
                           detector_star_pos_dither_pos_1["sigma"],
                           true_star_pos["epsilon"],
                           true_star_pos["sigma"])
    }


dither_pos_2_baselines = {} # none!


dither_pos_3_baselines = {
    
    "dp_3_del_zeta_eta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_3["zeta"], 
                           detector_star_pos_dither_pos_3["eta"],
                           true_star_pos["zeta"],
                           true_star_pos["eta"]),
    
    "dp_3_del_zeta_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_3["zeta"], 
                           detector_star_pos_dither_pos_3["theta"],
                           true_star_pos["zeta"],
                           true_star_pos["theta"]),

    "dp_3_del_zeta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_3["zeta"], 
                           detector_star_pos_dither_pos_3["sigma"],
                           true_star_pos["zeta"],
                           true_star_pos["sigma"]),
    
#--------------------------------------------------------
    "dp_3_del_eta_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_3["eta"], 
                           detector_star_pos_dither_pos_3["theta"],
                           true_star_pos["eta"],
                           true_star_pos["theta"]),

    "dp_3_del_eta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_3["eta"], 
                           detector_star_pos_dither_pos_3["sigma"],
                           true_star_pos["eta"],
                           true_star_pos["sigma"]),

#--------------------------------------------------------
    "dp_3_del_sigma_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_3["sigma"], 
                           detector_star_pos_dither_pos_3["theta"],
                           true_star_pos["sigma"],
                           true_star_pos["theta"])
    }


dither_pos_4_baselines = {
    
    "dp_4_del_D1_zeta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_4["D1"], 
                           detector_star_pos_dither_pos_4["zeta"],
                           true_star_pos["D1"],
                           true_star_pos["zeta"]),
    
    "dp_4_del_D1_eta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_4["D1"], 
                           detector_star_pos_dither_pos_4["eta"],
                           true_star_pos["D1"],
                           true_star_pos["eta"]),
     
    "dp_4_del_D1_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_4["D1"], 
                           detector_star_pos_dither_pos_4["sigma"],
                           true_star_pos["D1"],
                           true_star_pos["sigma"]),
    
#--------------------------------------------------------    
    "dp_4_del_zeta_eta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_4["zeta"], 
                           detector_star_pos_dither_pos_4["eta"],
                           true_star_pos["zeta"],
                           true_star_pos["eta"]),
    
    "dp_4_del_zeta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_4["zeta"], 
                           detector_star_pos_dither_pos_4["sigma"],
                           true_star_pos["zeta"],
                           true_star_pos["sigma"]),
    
#--------------------------------------------------------    
    "dp_4_del_eta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_4["eta"], 
                           detector_star_pos_dither_pos_4["sigma"],
                           true_star_pos["eta"],
                           true_star_pos["sigma"])
    }


dither_pos_5_baselines = {
    
    "dp_5_del_D1_zeta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_5["D1"], 
                           detector_star_pos_dither_pos_5["zeta"],
                           true_star_pos["D1"],
                           true_star_pos["zeta"]),
 
    "dp_5_del_D1_eta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_5["D1"], 
                           detector_star_pos_dither_pos_5["eta"],
                           true_star_pos["D1"],
                           true_star_pos["eta"]),
 
    "dp_5_del_D1_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_5["D1"], 
                           detector_star_pos_dither_pos_5["sigma"],
                           true_star_pos["D1"],
                           true_star_pos["sigma"]),
 
#--------------------------------------------------------    
    "dp_5_del_zeta_eta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_5["zeta"], 
                           detector_star_pos_dither_pos_5["eta"],
                           true_star_pos["zeta"],
                           true_star_pos["eta"]),
 
    "dp_5_del_zeta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_5["zeta"], 
                           detector_star_pos_dither_pos_5["sigma"],
                           true_star_pos["zeta"],
                           true_star_pos["sigma"]),
 
#-------------------------------------------------------- 
    "dp_5_del_eta_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_5["eta"], 
                           detector_star_pos_dither_pos_5["sigma"],
                           true_star_pos["eta"],
                           true_star_pos["sigma"])
    }



dither_pos_6_baselines = {
    
    "dp_6_del_D1_sigma": 
    plate_scale_and_angles(detector_star_pos_dither_pos_6["D1"], 
                           detector_star_pos_dither_pos_6["sigma"],
                           true_star_pos["D1"],
                           true_star_pos["sigma"])
    }


dither_pos_7_baselines = {

    "dp_7_del_D1_zeta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_7["D1"], 
                           detector_star_pos_dither_pos_7["zeta"],
                           true_star_pos["D1"],
                           true_star_pos["zeta"]),
  
    "dp_7_del_D1_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_7["D1"], 
                           detector_star_pos_dither_pos_7["theta"],
                           true_star_pos["D1"],
                           true_star_pos["theta"]),
  
#--------------------------------------------------------
    "dp_7_del_zeta_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_7["zeta"], 
                           detector_star_pos_dither_pos_7["theta"],
                           true_star_pos["zeta"],
                           true_star_pos["theta"])
    }




dither_pos_8_baselines = {} # none!


dither_pos_9_baselines = {
    
    "dp_9_del_zeta_eta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_9["zeta"], 
                           detector_star_pos_dither_pos_9["eta"],
                           true_star_pos["zeta"],
                           true_star_pos["eta"]),
  
    "dp_9_del_zeta_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_9["zeta"], 
                           detector_star_pos_dither_pos_9["theta"],
                           true_star_pos["zeta"],
                           true_star_pos["theta"]),
  
#--------------------------------------------------------
    "dp_9_del_eta_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_9["eta"], 
                           detector_star_pos_dither_pos_9["theta"],
                           true_star_pos["eta"],
                           true_star_pos["theta"])
    }


dither_pos_10_baselines = {
    
    "dp_10_del_B1_zeta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["B1"], 
                           detector_star_pos_dither_pos_10["zeta"],
                           true_star_pos["B1"],
                           true_star_pos["zeta"]),
  
    "dp_10_del_B1_eta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["B1"], 
                           detector_star_pos_dither_pos_10["eta"],
                           true_star_pos["B1"],
                           true_star_pos["eta"]),
  
    "dp_10_del_B1_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["B1"], 
                           detector_star_pos_dither_pos_10["theta"],
                           true_star_pos["B1"],
                           true_star_pos["theta"]),
  
    "dp_10_del_B1_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["B1"], 
                           detector_star_pos_dither_pos_10["A1"],
                           true_star_pos["B1"],
                           true_star_pos["A1"]),
  
    "dp_10_del_B1_B2": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["B1"], 
                           detector_star_pos_dither_pos_10["B2"],
                           true_star_pos["B1"],
                           true_star_pos["B2"]),
  
#--------------------------------------------------------      
    "dp_10_del_B2_zeta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["B2"], 
                           detector_star_pos_dither_pos_10["zeta"],
                           true_star_pos["B2"],
                           true_star_pos["zeta"]),
  
    "dp_10_del_B2_eta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["B2"], 
                           detector_star_pos_dither_pos_10["eta"],
                           true_star_pos["B2"],
                           true_star_pos["eta"]),
  
    "dp_10_del_B2_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["B2"], 
                           detector_star_pos_dither_pos_10["theta"],
                           true_star_pos["B2"],
                           true_star_pos["theta"]),
  
    "dp_10_del_B2_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["B2"], 
                           detector_star_pos_dither_pos_10["A1"],
                           true_star_pos["B2"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------         
    "dp_10_del_zeta_eta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["zeta"], 
                           detector_star_pos_dither_pos_10["eta"],
                           true_star_pos["zeta"],
                           true_star_pos["eta"]),
  
    "dp_10_del_zeta_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["zeta"], 
                           detector_star_pos_dither_pos_10["theta"],
                           true_star_pos["zeta"],
                           true_star_pos["theta"]),
  
    "dp_10_del_zeta_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["zeta"], 
                           detector_star_pos_dither_pos_10["A1"],
                           true_star_pos["zeta"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------    
    "dp_10_del_eta_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["eta"], 
                           detector_star_pos_dither_pos_10["theta"],
                           true_star_pos["eta"],
                           true_star_pos["theta"]),
  
    "dp_10_del_eta_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["eta"], 
                           detector_star_pos_dither_pos_10["A1"],
                           true_star_pos["eta"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------
    "dp_10_del_theta_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_10["theta"], 
                           detector_star_pos_dither_pos_10["A1"],
                           true_star_pos["theta"],
                           true_star_pos["A1"])
    }



dither_pos_11_baselines = {
    
    "dp_11_del_B1_B2": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["B1"], 
                           detector_star_pos_dither_pos_11["B2"],
                           true_star_pos["B1"],
                           true_star_pos["B2"]),
  
    "dp_11_del_B1_beta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["B1"], 
                           detector_star_pos_dither_pos_11["beta"],
                           true_star_pos["B1"],
                           true_star_pos["beta"]),
  
    "dp_11_del_B1_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["B1"], 
                           detector_star_pos_dither_pos_11["theta"],
                           true_star_pos["B1"],
                           true_star_pos["theta"]),
  
    "dp_11_del_B1_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["B1"], 
                           detector_star_pos_dither_pos_11["E1"],
                           true_star_pos["B1"],
                           true_star_pos["E1"]),
  
    "dp_11_del_B1_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["B1"], 
                           detector_star_pos_dither_pos_11["A1"],
                           true_star_pos["B1"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------         
    "dp_11_del_B2_beta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["B2"], 
                           detector_star_pos_dither_pos_11["beta"],
                           true_star_pos["B2"],
                           true_star_pos["beta"]),
  
    "dp_11_del_B2_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["B2"], 
                           detector_star_pos_dither_pos_11["theta"],
                           true_star_pos["B2"],
                           true_star_pos["theta"]),
  
    "dp_11_del_B2_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["B2"], 
                           detector_star_pos_dither_pos_11["E1"],
                           true_star_pos["B2"],
                           true_star_pos["E1"]),
  
    "dp_11_del_B2_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["B2"], 
                           detector_star_pos_dither_pos_11["A1"],
                           true_star_pos["B2"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------    
    "dp_11_del_beta_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["beta"], 
                           detector_star_pos_dither_pos_11["theta"],
                           true_star_pos["beta"],
                           true_star_pos["theta"]),
  
    "dp_11_del_beta_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["beta"], 
                           detector_star_pos_dither_pos_11["E1"],
                           true_star_pos["beta"],
                           true_star_pos["E1"]),
  
    "dp_11_del_beta_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["beta"], 
                           detector_star_pos_dither_pos_11["A1"],
                           true_star_pos["beta"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------
    "dp_11_del_theta_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["theta"], 
                           detector_star_pos_dither_pos_11["E1"],
                           true_star_pos["theta"],
                           true_star_pos["E1"]),
  
    "dp_11_del_theta_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["theta"], 
                           detector_star_pos_dither_pos_11["A1"],
                           true_star_pos["theta"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------
    "dp_11_del_E1_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_11["E1"], 
                           detector_star_pos_dither_pos_11["A1"],
                           true_star_pos["E1"],
                           true_star_pos["A1"])
    }



dither_pos_12_baselines = {
    
    "dp_12_del_B1_B2": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["B1"], 
                           detector_star_pos_dither_pos_12["B2"],
                           true_star_pos["B1"],
                           true_star_pos["B2"]),
  
    "dp_12_del_B1_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["B1"], 
                           detector_star_pos_dither_pos_12["theta"],
                           true_star_pos["B1"],
                           true_star_pos["theta"]),
  
    "dp_12_del_B1_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["B1"], 
                           detector_star_pos_dither_pos_12["E1"],
                           true_star_pos["B1"],
                           true_star_pos["E1"]),
  
    "dp_12_del_B1_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["B1"], 
                           detector_star_pos_dither_pos_12["A1"],
                           true_star_pos["B1"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------    
    "dp_12_del_B2_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["B2"], 
                           detector_star_pos_dither_pos_12["theta"],
                           true_star_pos["B2"],
                           true_star_pos["theta"]),
  
    "dp_12_del_B2_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["B2"], 
                           detector_star_pos_dither_pos_12["E1"],
                           true_star_pos["B2"],
                           true_star_pos["E1"]),
  
    "dp_12_del_B2_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["B2"], 
                           detector_star_pos_dither_pos_12["A1"],
                           true_star_pos["B2"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------
    "dp_12_del_theta_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["theta"], 
                           detector_star_pos_dither_pos_12["E1"],
                           true_star_pos["theta"],
                           true_star_pos["E1"]),
  
    "dp_12_del_theta_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["theta"], 
                           detector_star_pos_dither_pos_12["A1"],
                           true_star_pos["theta"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------
    "dp_12_del_E1_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_12["E1"], 
                           detector_star_pos_dither_pos_12["A1"],
                           true_star_pos["E1"],
                           true_star_pos["A1"])
    }



dither_pos_13_baselines = {
    
    "dp_13_del_B1_B2": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["B1"], 
                           detector_star_pos_dither_pos_13["B2"],
                           true_star_pos["B1"],
                           true_star_pos["B2"]),
  
    "dp_13_del_B1_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["B1"], 
                           detector_star_pos_dither_pos_13["theta"],
                           true_star_pos["B1"],
                           true_star_pos["theta"]),
  
    "dp_13_del_B1_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["B1"], 
                           detector_star_pos_dither_pos_13["E1"],
                           true_star_pos["B1"],
                           true_star_pos["E1"]),
  
    "dp_13_del_B1_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["B1"], 
                           detector_star_pos_dither_pos_13["A1"],
                           true_star_pos["B1"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------    
    "dp_13_del_B2_theta": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["B2"], 
                           detector_star_pos_dither_pos_13["theta"],
                           true_star_pos["B2"],
                           true_star_pos["theta"]),
  
    "dp_13_del_B2_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["B2"], 
                           detector_star_pos_dither_pos_13["E1"],
                           true_star_pos["B2"],
                           true_star_pos["E1"]),
  
    "dp_13_del_B2_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["B2"], 
                           detector_star_pos_dither_pos_13["A1"],
                           true_star_pos["B2"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------
    "dp_13_del_theta_E1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["theta"], 
                           detector_star_pos_dither_pos_13["E1"],
                           true_star_pos["theta"],
                           true_star_pos["E1"]),
  
    "dp_13_del_theta_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["theta"], 
                           detector_star_pos_dither_pos_13["A1"],
                           true_star_pos["theta"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------
    "dp_13_del_E1_A1": 
    plate_scale_and_angles(detector_star_pos_dither_pos_13["E1"], 
                           detector_star_pos_dither_pos_13["A1"],
                           true_star_pos["E1"],
                           true_star_pos["A1"])
    }



dither_pos_14_baselines = {

    "dp_14_del_B1_B2": 
    plate_scale_and_angles(detector_star_pos_dither_pos_14["B1"], 
                           detector_star_pos_dither_pos_14["B2"],
                           true_star_pos["B1"],
                           true_star_pos["B2"]),
  
    "dp_14_del_B1_E1":
    plate_scale_and_angles(detector_star_pos_dither_pos_14["B1"], 
                           detector_star_pos_dither_pos_14["E1"],
                           true_star_pos["B1"],
                           true_star_pos["E1"]),
  
    "dp_14_del_B1_A1":
    plate_scale_and_angles(detector_star_pos_dither_pos_14["B1"], 
                           detector_star_pos_dither_pos_14["A1"],
                           true_star_pos["B1"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------
    "dp_14_del_B2_E1":
    plate_scale_and_angles(detector_star_pos_dither_pos_14["B2"], 
                           detector_star_pos_dither_pos_14["E1"],
                           true_star_pos["B2"],
                           true_star_pos["E1"]),
  
    "dp_14_del_B2_A1":
    plate_scale_and_angles(detector_star_pos_dither_pos_14["B2"], 
                           detector_star_pos_dither_pos_14["A1"],
                           true_star_pos["B2"],
                           true_star_pos["A1"]),
  
#--------------------------------------------------------
    "dp_14_del_E1_A1":
    plate_scale_and_angles(detector_star_pos_dither_pos_14["E1"], 
                           detector_star_pos_dither_pos_14["A1"],
                           true_star_pos["E1"],
                           true_star_pos["A1"])
    }



dither_pos_15_baselines = {

    "dp_15_del_B1_B2":
    plate_scale_and_angles(detector_star_pos_dither_pos_15["B1"], 
                           detector_star_pos_dither_pos_15["B2"],
                           true_star_pos["B1"],
                           true_star_pos["B2"]),
  
    "dp_15_del_B1_E1":
    plate_scale_and_angles(detector_star_pos_dither_pos_15["B1"], 
                           detector_star_pos_dither_pos_15["E1"],
                           true_star_pos["B1"],
                           true_star_pos["E1"]),
  
#--------------------------------------------------------
    "dp_15_del_B2_E1":
    plate_scale_and_angles(detector_star_pos_dither_pos_15["B2"], 
                           detector_star_pos_dither_pos_15["E1"],
                           true_star_pos["B2"],
                           true_star_pos["E1"])
    }



dither_pos_16_baselines = {

    "dp_16_del_B2_E1":
    plate_scale_and_angles(detector_star_pos_dither_pos_16["B2"], 
                           detector_star_pos_dither_pos_16["E1"],
                           true_star_pos["B2"],
                           true_star_pos["E1"])
    }



dither_pos_17_baselines = {} # none!



# print stuff FYI
for item in dither_pos_0_baselines:
    print(item, dither_pos_0_baselines[item])
for item in dither_pos_1_baselines:
    print(item, dither_pos_1_baselines[item])
for item in dither_pos_2_baselines:
    print(item, dither_pos_2_baselines[item])
for item in dither_pos_3_baselines:
    print(item, dither_pos_3_baselines[item])
for item in dither_pos_4_baselines:
    print(item, dither_pos_4_baselines[item])
for item in dither_pos_5_baselines:
    print(item, dither_pos_5_baselines[item])
for item in dither_pos_6_baselines:
    print(item, dither_pos_6_baselines[item])
for item in dither_pos_7_baselines:
    print(item, dither_pos_7_baselines[item])
for item in dither_pos_8_baselines:
    print(item, dither_pos_8_baselines[item])
for item in dither_pos_9_baselines:
    print(item, dither_pos_9_baselines[item])
for item in dither_pos_10_baselines:
    print(item, dither_pos_10_baselines[item])
for item in dither_pos_11_baselines:
    print(item, dither_pos_11_baselines[item])
for item in dither_pos_12_baselines:
    print(item, dither_pos_12_baselines[item])
for item in dither_pos_13_baselines:
    print(item, dither_pos_13_baselines[item])
for item in dither_pos_14_baselines:
    print(item, dither_pos_14_baselines[item])
for item in dither_pos_15_baselines:
    print(item, dither_pos_15_baselines[item])
for item in dither_pos_16_baselines:
    print(item, dither_pos_16_baselines[item])
for item in dither_pos_17_baselines:
    print(item, dither_pos_17_baselines[item])


########################################################
# PUT PLATE SCALE, ANGLE VALUES INTO CONVENIENT ARRAYS

baselineNameArray = []
plateScaleArray = []
angularChangeArray = []
plateScaleElem = 0 # index for the element denoting angle differences
angleElem = 1

for key in dither_pos_0_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_0_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_0_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_0_baselines[item][angleElem]])

for key in dither_pos_1_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)    
for item in dither_pos_1_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_1_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_1_baselines[item][angleElem]])

for key in dither_pos_2_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_2_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_2_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_2_baselines[item][angleElem]])

for key in dither_pos_3_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_3_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_3_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_3_baselines[item][angleElem]])

for key in dither_pos_4_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)    
for item in dither_pos_4_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_4_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_4_baselines[item][angleElem]])

for key in dither_pos_5_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_5_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_5_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_5_baselines[item][angleElem]])

for key in dither_pos_6_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_6_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_6_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_6_baselines[item][angleElem]])

for key in dither_pos_7_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_7_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_7_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_7_baselines[item][angleElem]])

for key in dither_pos_8_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_8_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_8_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_8_baselines[item][angleElem]])

for key in dither_pos_9_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_9_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_9_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_9_baselines[item][angleElem]])

for key in dither_pos_10_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_10_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_10_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_10_baselines[item][angleElem]])

for key in dither_pos_11_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_11_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_11_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_11_baselines[item][angleElem]])

for key in dither_pos_12_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_12_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_12_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_12_baselines[item][angleElem]])

for key in dither_pos_13_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_13_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_13_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_13_baselines[item][angleElem]])

for key in dither_pos_14_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_14_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_14_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_14_baselines[item][angleElem]])

for key in dither_pos_15_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_15_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_15_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_15_baselines[item][angleElem]])

for key in dither_pos_16_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_16_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_16_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_16_baselines[item][angleElem]])

for key in dither_pos_17_baselines:
    baselineNameArray = np.append([baselineNameArray],
                                  key)
for item in dither_pos_17_baselines:
    plateScaleArray = np.append([plateScaleArray],
                          [dither_pos_17_baselines[item][plateScaleElem]])
    angularChangeArray = np.append([angularChangeArray],
                          [dither_pos_17_baselines[item][angleElem]])

# issue with some stars in pairs being of almost same DEC, where code thinks one star is further N
# than the other when in 'true' coordinates, but the other star in LMIRCam pixel coordinates
if (len(np.where(angularChangeArray < -170)) > 0): # if there are angle differences close to -180
    print('---------------------------------------------------------------------------')
    print('The following angle differences are being corrected for a 180 deg flip...')
    print('(LMIRCam thinks the rightmost star is higher, true coordinates say the opposite)')
    print(angularChangeArray[np.where(angularChangeArray < -170)])
    angularChangeArray[np.where(angularChangeArray < -170)] += 180

# ... and if the reverse is the case
if (len(np.where(angularChangeArray < -170)) < 0): # if there are angle differences close to -180
    print('---------------------------------------------------------------------------')
    print('The following angle differences are being corrected for a 180 deg flip...')
    print('(LMIRCam thinks the rightmost star is lower, true coordinates say the opposite)')
    print(angularChangeArray[np.where(angularChangeArray > 170)])
    angularChangeArray[np.where(angularChangeArray > 170)] -= 180

########################################################
# OPTIONAL: INVESTIGATE TRENDS IN THE ARRAYS (use manually)

# indices that sort by angular change
#indicesSorted = np.argsort(angularChangeArray)
#baselineNamesSorted = baselineNameArray[indicesSorted]
#plateScaleSorted = plateScaleArray[indicesSorted]
#angularChangeSorted = angularChangeArray[indicesSorted]

    
########################################################
# CALCULATE FINAL RESULTS

# decimals to represent 1-sigma bounds
sigmaBoundLow = 0.15865
sigmaBoundHigh = 0.84135

# angle offset
medianAngle = np.percentile(angularChangeArray, 50, interpolation='linear')
lowerEndAngle = np.subtract(np.percentile(angularChangeArray, 50, interpolation='linear'),
                            np.percentile(angularChangeArray, 100.*sigmaBoundLow, interpolation='linear'))
upperEndAngle = np.subtract(np.percentile(angularChangeArray, 100.*sigmaBoundHigh, interpolation='linear'),
                            np.percentile(angularChangeArray, 50, interpolation='linear'))

# plate scale
medianPlateScale = 1000.*np.percentile(plateScaleArray, 50, interpolation='linear')
lowerEndPlateScale = 1000.*np.subtract(np.percentile(plateScaleArray, 50, interpolation='linear'),
                                       np.percentile(plateScaleArray, 100.*sigmaBoundLow, interpolation='linear'))
upperEndPlateScale = 1000.*np.subtract(np.percentile(plateScaleArray, 100.*sigmaBoundHigh, interpolation='linear'),
                                       np.percentile(plateScaleArray, 50, interpolation='linear'))

# print results

print('----------------------')
print('Angular offset E of N:')
print(medianAngle)
print('+'+str(upperEndAngle))
print('-'+str(lowerEndAngle))

print('----------------------')
print('Plate scale (mas/pix):')
print(medianPlateScale)
print('+'+str(upperEndPlateScale))
print('-'+str(lowerEndPlateScale))


# plot

plt.figure(figsize=(16,8))

plt.subplot(121)
plt.hist(1000.*plateScaleArray, bins=25, color='grey')
plt.axvline(medianPlateScale-lowerEndPlateScale, linewidth=5, linestyle='dashed', color='k', snap=False)
plt.axvline(medianPlateScale+upperEndPlateScale, linewidth=5, linestyle='dashed', color='k', snap=False)
plt.title("LMIRCam Plate Scale from Star Baseline Separations\n(from Trapezium data taken UT 2016 Nov 12)")
plt.xlabel("Plate Scale (mas/pix)")
plt.ylabel("Number of baselines in bin")

plt.subplot(122)
plt.hist(angularChangeArray, bins=25, color='grey')
plt.axvline(medianAngle-lowerEndAngle, linewidth=5, linestyle='dashed', color='k', snap=False)
plt.axvline(medianAngle+upperEndAngle, linewidth=5, linestyle='dashed', color='k', snap=False)
plt.title("LMIRCam Orientation from Star Baseline Angles\n(from Trapezium data taken UT 2016 Nov 12)")
plt.xlabel("Detector vertical axis angle E of true N (deg)")
plt.ylabel("Number of baselines in bin")

plt.show()
