Introduction
=================
This software is for dewarping the illumination beam on a detector
based on a known 'perfect' distribution of point source, and also for
finding the rotation solution of the dewarped frames based on the
known positions of stellar objects.

This was originally made for the LBTI LMIRcam detector, but is
generalizeable to any project with frames with a grid of point
sources, an astrometric field.

The pinhole images are used to make the dewarping solution.

The astrometric field is used to determine the on-sky orientation and
plate scale. (As of now, it does not play a role in the dewarping solution.)

We want to find the polynomial coefficients that map between empirical
pinhole locations and an idealized grid, namely

(eqns. 1 and 2 in procedure)
