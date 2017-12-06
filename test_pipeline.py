# This is the script that should run the entire pipeline to reduce the delta Cyg data from Feb 2017

# created by E.S., 2017 Nov 11

from astrom import match_pinholes
from astrom import find_dewarp_solution, apply_dewarp_solution

filenameString = 'test2' # this string is added onto pickle files
match_pinholes.match_pinholes(
    [110,-110],
    48.0,
    [100,512],
    -1e-9,
    -3.2,
    filenameString,
    plot=False) # do the reduction

find_dewarp_solution.find_dewarp(
    filenameString)

apply_dewarp_solution.apply_dewarp(8597,8942,filenameString)

