Dewarping
=================
This software is for dewarping the illumination beam on a detector
based on a known 'perfect' distribution of point source, and also for
finding the rotation solution of the dewarped frames based on the
known positions of stellar objects.

make dewarp coords()

Within the script find dewarp solution.py, you will see some functions and arrays appear in the first section of the code that you may have to run through a couple times so that you can tweak the function inputs to values that are optimal for your data. (See comments in the code for details.) It’s also prob- ably good to mask pinholes in the heavily vignetted region of the array (Fig. 2).

Once that’s done, run the script again so that it runs past the
function make dewarp coords(). This finds the aforementioned
coefficients by solving a least-squares problem via Moore-Penrose
pseudoinverse matrices. (I find J. Stone’s condensed description of
this to be helpful.) Schematically, what is being done is shown in
cartoon form in Fig. 3.

(Fig. 3 in procedure)
.. image:: images/kinks2.png
	   :width: 400
	   :alt: Alternative text

dewarp with precomputed coords()

This next function takes the raw image, pastes the warped coordinates onto it, and then smooths everything out by resampling the image point-by-point over the entire image space, interpolating as needed when the coordinates are not at integer values (Fig. 4).
As a check, closely compare the pinhole grid images before and after (Fig. 5).

(Fig. 4 in procedure)

The last part of the script makes a lovely barb plot, putting evenly-spaced vectors over the array to show the directions that points on the readouts have to be stretched in order to dewarp it (Fig. 6).