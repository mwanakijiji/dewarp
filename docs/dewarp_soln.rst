Dewarping
=================

**********
Requirements
**********

A FITS frame of sources with patterned, known locations (like
a pinhole grid in Fig. :numref:`pinhole_ex`).


**********
The idea
**********

This software is for dewarping the illumination beam on a detector
based on a known 'perfect' distribution of point source, and also for
finding the rotation solution of the dewarped frames based on the
known positions of stellar objects.

:math:`x_{i}=\sum^{N}_{i=0}\sum^{N}_{j=0}K_{x}^{(i,j)}x_{o}^{(j)}y_{o}^{(i)}`

:math:`\underbrace{y_{i}}_\text{warped}=\sum^{N}_{i=0}\sum^{N}_{j=0}K_{y}^{(i,j)}\underbrace{x_{o}^{(j)}y_{o}^{(i)}}_\text{dewarped}`


$K_{x}^{(i,j)}$ and $K_{y}^{(i,j)}$
      
make dewarp coords()

Within the script find dewarp solution.py, you will see some functions and arrays appear in the first section of the code that you may have to run through a couple times so that you can tweak the function inputs to values that are optimal for your data. (See comments in the code for details.) It’s also prob- ably good to mask pinholes in the heavily vignetted region of the array (Fig. 2).

Once that’s done, run the script again so that it runs past the
function make dewarp coords(). This finds the aforementioned
coefficients by solving a least-squares problem via Moore-Penrose
pseudoinverse matrices. (I find J. Stone’s condensed description of
this to be helpful.) Schematically, what is being done is shown in
cartoon form in Fig. 3.

(Fig. 3 in procedure)

.. _label:
.. figure:: images/kinks2.png
	   :scale: 90 %
           :align: center
	   :alt: Alternative text

dewarp with precomputed coords()

This next function takes the raw image, pastes the warped coordinates onto it, and then smooths everything out by resampling the image point-by-point over the entire image space, interpolating as needed when the coordinates are not at integer values (Fig. 4).
As a check, closely compare the pinhole grid images before and after (Fig. 5).

(Fig. 4 in procedure)

The last part of the script makes a lovely barb plot, putting evenly-spaced vectors over the array to show the directions that points on the readouts have to be stretched in order to dewarp it (Fig. 6).
