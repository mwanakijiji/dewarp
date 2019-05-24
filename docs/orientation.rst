Sky orientation and plate scale
=================

**********
Requirements
**********

#. Stellar astrometric images (preferably dewarped)
#. A csv of known stellar targets with RA and DEC

**********
The idea
**********

We find star positions in pixel space, and compare their locations in
RA, DEC. All possible baselines are used between identified stars to
find their separations and angles relative to north at PA=0.

**********
Procedure
**********

(Fig. 6 in pdf)

previously-constrained ‘true’ astrometric values to find the plate
scale, and the angles are compared with their ‘true’ counterparts to
find the residual angular offset of the detector. Note that in one
frame of Ns stars, the total number of baselines among stellar pairs
is “Ns choose 2”:

:math:`N_{b}=\begin{pmatrix}N_{s} \\ 2\end{pmatrix} \equiv \frac{N_{s}!}{2!(N_{s}-2)!}=\frac{N_{s}(N_{s}-1)}{2}`

2.2 Find stars in pixel space: find asterism star locations.py

The pre-requisite dewarping of images is performed in find dewarp
solution.py. But you also need to derotate them based on their
parallactic angle. For this I used a quick IDL script that read the
headers and used the IDL function rot. (I included a FYI copy of my
script derotate trapezium data ut 2016 11 12.pro.) However you do it,
take the dewarped FITS files (which are now residing in /step02
dewarped/, derotate them, and write out the results to the
directory /step03 derotate/.

Now take the median of each dither position. Dump these medians into /step04 ditherMedians/. (Since calibration-related images are often taken during mediocre conditions, it
is important to take many images at a given dither position.) Overlay the resul-
tant dither medians (Adobe Photoshop is very useful for this), mark any visible

stars, and then cross-check them with a known astrometric source. In
the case of the Trapezium Cluster, one can use the images and Table 1
in Close+ 2012 ApJ 749:180. In Fig. 7, I labeled stars using the
conventions in Close+ 2012, and used my own Greek lettering if they
were without label.

(Fig. 7 in procedure)

The script find asterism star locations.py takes the intermediary step of
determining star locations in pixel space, and printing locations in pixel space to
the screen. Check each centroid manually in the plot, to see if it’s a real star or
not (Fig. 8). Copy the true positive locations in pixel space that are returned in
the Terminal, and populate the dictionaries in the script find plate
scale and orientation.py.

Find plate scale, angle offset: find plate scale and orientation.py

(Fig. 12 in pdf)
