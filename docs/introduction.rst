Introduction
=================
This software is for

#. determining and removing detector distortion
#. finding a correct detector orientation relative to true North
#. determining the detector plate scale.

This was originally made for the LBTI LMIRcam detector, but is
generalizeable to any project with reference frames of ideal sources
and/or an astrometric field paired with a csv of astrometric target locations.
   
The distortion removal is based on a mapping made from a distribution
of empirical sources (like a pinhole grid) and an idealized
distribution of point sources with a "perfect" distortion. This is
them mapped to a perfect distribition without distortions.

The detector orientation and plate scale is determined by allowing the user to mark
the identities of stars in a series of astrometric field
frames. Following this, all lengths and angles of non-redundant baselines between pairs of
stars is calculated to find both the plate scale and angle relative to
true North.

Attribution
-----------

Please cite `Foreman-Mackey, Hogg, Lang & Goodman (2012)
<http://arxiv.org/abs/1202.3665>`_ if you find this code useful in your
research and add your paper to `the testimonials list
<https://github.com/dfm/emcee/blob/master/docs/testimonials.rst>`_.
The BibTeX entry for the paper is::

    @article{emcee,
       author = {{Foreman-Mackey}, D. and {Hogg}, D.~W. and {Lang}, D. and {Goodman}, J.},
        title = {emcee: The MCMC Hammer},
      journal = {PASP},
         year = 2013,
       volume = 125,
        pages = {306-312},
       eprint = {1202.3665},
          doi = {10.1086/670067}
    }
