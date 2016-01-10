# wht_reduction_scripts
Reductions scripts for WHT using specreduce

## Installation Instructions

The following scripts require these packages to be installed in order to run.  These can all be installed via pip:

+ numpy
+ scipy
+ matplotlib
+ astropy
+ astroscrappy
+ pyspectrorgaph


At this time, the development versions of these two packages need to be installed
+ [ccdproc](https://github.com/astropy/ccdproc.git)
+ [specreduce](https://github.com/crawfordsm/specreduce.git)

The specreduce does depend on [PyQt4](https://riverbankcomputing.com/software/pyqt/intro) package. 

## Instructions

To perform basic data reductions, follow these steps:


1. Run the basic reductions on the data.   Provide the full path to the directory with the raw data and to the directory where the product data should be produced.  These should not be the same directory.

    python wht_basic_reductions.py [full_path_to_raw_data] [full_path_to_reduced_data]

2. Create wave maps for each of the arc frames

3.  Apply the wavemaps to object frames


