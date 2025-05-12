## Preparing data

The codes in this directory prepare 5-min omni data for use in the TA16 model. To use, compile the *calculate_tilt.f* file using your favorite FORTRAN compiler and edit *calculate_tile.f*, *fill_IMF_and_SW_gaps.py*, and *calculate_TA16_inds.py* to point at the file of interest. Then simply run then in order,

 -1 *Fill_IMF_and_SW_gaps.py*
 -2 *calculate_tilt.out* (or whatever the compiled file is called)
 -3 *calculate_TA16_inds.py*

The formats of each of these files are contained in the *formats/* directory.
