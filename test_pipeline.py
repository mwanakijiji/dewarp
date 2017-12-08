# This is the script that should run the entire pipeline to reduce the delta Cyg data from Feb 2017

# created by E.S., 2017 Nov 11

from astrom import match_pinholes
from astrom import find_dewarp_solution, apply_dewarp_solution

filenameString = '20171108_both_apertures' # this string is added onto pickle files
plotTitleString = '2017 Nov 08 (both apertures)' # for distortion vector plot title

# match the empirical and ideal pinholes
match_pinholes.match_pinholes(
    [110,-110],
    48.0,
    [100,512],
    -1e-9,
    -3.2,
    filenameString,
    plotTitleString,
    plot=True) # add another flag for when user is just testing model grid first

# work out the solution
find_dewarp_solution.find_dewarp(
    filenameString,
    plotTitleString,
    plot=True)

# apply solution
apply_dewarp_solution.apply_dewarp(
    8597,
    8942,
    filenameString)

