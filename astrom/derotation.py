import numpy as np
import configparser
import glob
from astropy.io import fits
from astrom_lmircam_soln import *
from pytictoc import TicToc

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
config.read("astrom/config.ini")

def derotate_image_forloop(dateString):

    # obtain the list of files which have been dewarped, and need to de-rotated
    # according to the PA in the FITS headers
    asterism_frames_directory_retrieve = str(config["data_dirs"]["DIR_ASTERISM_DEWARP"])
    asterism_frames_pre_derot_names = list(glob.glob(os.path.join(asterism_frames_directory_retrieve, "*.fits")))

    for f in range(0,len(asterism_frames_pre_derot_names)): # loop over filenames
        print('----------------------------------------')

        t = TicToc() # create instance of timer
        t.tic() # start timer

        # retrieve image
        image, header = fits.getdata(asterism_frames_pre_derot_names[f],0,header=True)
        #print(header)

        try:
            # find PA from header
            pa = np.float(header['LBT_PARA'])
    
            # derotate
            image_derot = rot(im = image,
                              angle = -pa,
                              axis = [1024,1024],
                              order = 3, pivot=False) # axis coord here is just a dummy

            #print(angle)
            print('Derotating image ' + str(os.path.basename(asterism_frames_pre_derot_names[f])) + '...')

            print('ya')
            # save
            fits.writeto(str(config["data_dirs"]["DIR_ASTERISM_DEROT"] + \
                         "derotated_" + \
                         os.path.basename(asterism_frames_pre_derot_names[f])),
                         image_derot, header, overwrite=False)

            t.toc()
            print('------------------------------')

        except:
            print("Frame " + str(os.path.basename(asterism_frames_pre_derot_names[f])) + \
                  " has no parallactic angle. Maybe header keyword being sought is wrong?")
