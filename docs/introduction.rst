Introduction
=================
This software is for

#. determining and removing detector distortion
#. finding a correct detector orientation relative to true North
#. determining the detector plate scale

This was originally made for the LBTI LMIRcam detector, but is
generalizable to any project with reference frames of ideal sources
and/or an astrometric field paired with a .csv file of astrometric target locations.

The distortion removal is based on a mapping made from a distribution
of empirical sources (like a pinhole grid) and an idealized
distribution of point sources with a "perfect" distortion. This is
them mapped to a perfect distribution without distortions.

The detector orientation and plate scale is determined by allowing the user to mark
the identities of stars in a series of astrometric field
frames. Following this, all lengths and angles of non-redundant baselines between pairs of
stars is calculated to find both the plate scale and angle relative to
true North.

Attribution
-----------

Please cite Dewarp `Spalding and Stone (2019)` if you used this code or distortion
corrections generated from it.
The BibTeX entry is::

    @misc{dewarp,
      author = {{Spalding}, E. and J. {Stone}},
      title = {Dewarp},
      keywords = {Software},
      howpublished = {Astrophysics Source Code Library},
      year = 2019,
      archivePrefix = "ascl",
      eprint = {1907.008},
      adsurl = {https://ui.adsabs.harvard.edu/abs/2019ascl.soft07008S/abstract},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
