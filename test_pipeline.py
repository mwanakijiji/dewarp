# This is the script that should run the entire pipeline to reduce the delta Cyg data from Feb 2017

# created by E.S., 2017 Nov 11

from astrom import match_pinholes
from astrom import find_dewarp_solution, apply_dewarp_solution, derotation, intermission, find_asterism_star_locations, comparison
import ipdb

dateStringShort = '190125'
filenameString = dateStringShort+'_DX_only' # this string is added onto pickle files
plotTitleString = 'Pinholes 2019 Jan 25 (DX only)' # for plot titles

# frame numbers of frames that will need to be derotated and dewarped
startNum = 8597
stopNum = 8941 # inclusive

### PART 1: FIND DEWARP SOLUTION
# match the empirical and ideal pinholes
# barrelCenterPass is (x,y)

match_pinholes.match_pinholes(
    translationPass = [0,-10],
    holeSpacingPass = 48.0,
    barrelCenterPass = [1124,715],
    barrelAmountPass = 0,
    rotationAnglePass = -0.8,
    writeoutString = filenameString,
    plotTitleString = plotTitleString,
    plot=True) # add another flag for when user is just testing model grid first

# work out the solution
find_dewarp_solution.find_dewarp(
    filenameString,
    plotTitleString,
    plot=True)

'''
### PART 2: APPLY DEWARP AND DEROTATE

# apply dewarp solution
apply_dewarp_solution.apply_dewarp(
    startNum,
    stopNum,
    filenameString,
    maskUnsampled=True)

# derotate
derotation.derotate_image_forloop(startNum,stopNum,dateStringShort)


### INTERMISSION: USER NEEDS TO MAKE DITHER MEDIANS
intermission.prompt_make_dithers()


### PART 3: FIND PLATE SCALE

# identify asterism stars in (x,y) space
find_asterism_star_locations.find_stars(
    dateStringShort,
    17)

# find position angle offset and plate scale
comparison.angOffset_plateScale(dateStringShort,
                                plotTitleString,
                                plot=True)
'''
