# This is the script that should run the entire pipeline to reduce the delta Cyg data from Feb 2017

# created by E.S., 2017 Nov 11

from astrom import match_pinholes
from astrom import find_dewarp_solution

match_pinholes.match_pinholes(
    translation=[110,-110],
    holeSpacing=48.0,
    barrelCenter=[100,512],
    barrelAmount=-1e-9,
    rotationAngle=-3.2,
    plot=False) # do the reduction

find_dewarp_solution()

