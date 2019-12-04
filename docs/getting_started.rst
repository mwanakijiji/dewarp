Getting started
=================
This repository provides importable functions that can be used to
assemble a pipeline as needed. An example is in the script
template_pipeline.py, which you can follow in the next
pages.

In the config.ini file, replace the DIR_HOME string with the path to your copy
of the repository. (The other directories will be made by the pipeline.)

Specify the the other directories in the config.ini file, such as that which
will contain the pinhole image in the config file.

Run

.. code-block:: bash

  python template_pipeline.py

which will proceed to make directories with the function make_dirs() and Then
get stuck, because it cannot find the FITS file of the pinholes. Put the pinhole
image into its directory, and rerun the above command.
