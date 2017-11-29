# created by E.S., 22 Nov 2016

import os

path=os.path.dirname(__file__)

dateStringAsterism = 'ut_2016_11_12'
dateStringPinholes = 'ut_2016_11_22'
generalStem = os.path.expanduser('~')+'/../../media/unasemaje/Seagate Expansion Drive/lbti_data_reduction/lmircam_astrometry/'

# data paths: trapezium
raw_trapezium_data_stem = generalStem+dateStringAsterism+'/asterism/rawData/' # obtain raw science data only (don't save anything to this!)
calibrated_trapezium_data_stem = generalStem+dateStringAsterism+'/asterism/processedData/' # deposit science arrays after bias-subtraction, flat-fielding, etc.
save_trapezium_data_stem = generalStem+'/textfile_results/' # deposit polynomial fit data

# data paths: pinholes
raw_pinholes_data_stem = generalStem+dateStringPinholes+'/pinholeGrid/rawData/' # obtain raw science data only (don't save anything to this!)
calibrated_pinholes_data_stem = generalStem+dateStringPinholes+'/pinholeGrid/processedData/' # deposit science arrays after bias-subtraction, flat-fielding, etc.
save_pinholes_data_stem = generalStem+'/textfile_results/' # deposit polynomial fit data
