## Preparing data

The codes in this directory prepare omni data for use in the TA16 model. To use, compile the *calculate_tilt.f* file using your favorite FORTRAN compiler and edit *fill_gaps.py* and *calculate_TA16_inds.py* to point at the files of interest. Then simply run each file in order,

 -1 *fill_gaps.py*
 -2 *calculate_tilt.out* \[list of dates as command line arg ex. 2024-01-23 for the file 2024-01-23_filled_gaps.list\] (or whatever the compiled file is called)
 -3 *calculate_TA16_inds.py*

The formats of each of these files are contained in the *formats/* directory.
