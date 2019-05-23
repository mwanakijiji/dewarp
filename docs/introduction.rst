Introduction
=================
This software is for
#1. removing detector distortion
#2. finding a correct detector orientation relative to true North
#3. determining the detector plate scale.

This was originally made for the LBTI LMIRcam detector, but is
generalizeable to any project with reference frames of ideal sources,
and/or an astrometric field paired with a csv of astrometric target locations.
   
The distortion removal is based on a mapping made from a distribution
of empirical sources (like a pinhole grid) and an idealized
distribution of point sources with a "perfect" distortion. This is
them mapped to a perfect distribition without distortions.

The detector orientation and plate scale is determined by allowing the user to mark
the identities of stars in a series of astrometric field
frames. Following this, all non-redundant baselines between pairs of
stars is calculated in both arcseconds on sky and in pixels on the detector.



The astrometric field is used to determine the on-sky orientation and
plate scale. 

We want to find the polynomial coefficients that map between empirical
pinhole locations and an idealized grid, namely

(eqns. 1 and 2 in procedure)
