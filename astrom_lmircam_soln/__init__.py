# created by E.S., 22 Nov 2016

import os
import subprocess

path=os.path.dirname(__file__)

dateStringAsterism = 'ut_2017_11_06'
dateStringPinholes = 'ut_2017_11_08'
generalStem = os.path.expanduser('~')+'/../../media/unasemaje/Seagate Expansion Drive/lbti_data_reduction/lmircam_astrometry/'

# data paths: trapezium
raw_trapezium_data_stem = generalStem+dateStringAsterism+'/asterism/rawData/' # obtain raw science data only (don't save anything to this!)
calibrated_trapezium_data_stem = generalStem+dateStringAsterism+'/asterism/processedData/' # deposit science arrays after bias-subtraction, flat-fielding, etc.
save_trapezium_data_stem = generalStem+'/textfile_results/' # deposit polynomial fit data

# data paths: pinholes
raw_pinholes_data_stem = generalStem+dateStringPinholes+'/pinholeGrid/rawData/' # obtain raw science data only (don't save anything to this!)
calibrated_pinholes_data_stem = generalStem+dateStringPinholes+'/pinholeGrid/processedData/step02_medianed/pinhole_image_median_bothApertures_171108.fits' # deposit science arrays after bias-subtraction, flat-fielding, etc.
save_pinholes_data_stem = generalStem+'/textfile_results/' # deposit polynomial fit data


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
