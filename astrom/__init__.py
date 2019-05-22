'''
Initialization
'''

import os
import numpy as np
import scipy
from scipy import ndimage, sqrt, stats, misc, signal
import git
import configparser
import multiprocessing


## SOME VARIABLES
# number of CPUs for parallelization
ncpu = multiprocessing.cpu_count()

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
config.read("astrom/config.ini")

def make_dirs():
    '''
    Make directories for housing files/info if they don't already exist
    '''

    # loop over all directory paths we will need
    for vals in config["data_dirs"]:
        abs_path_name = str(config["data_dirs"][vals])
        print("Directory exists: " + abs_path_name)

        # if directory does not exist, create it
        if not os.path.exists(abs_path_name):
            os.makedirs(abs_path_name)
            print("Made directory " + abs_path_name)
