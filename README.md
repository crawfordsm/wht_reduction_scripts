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

### Run the basic reductions on the data.   

Provide the full path to the directory with the raw data and to the directory where the product data should be produced.  These should not be the same directory.

    python wht_basic_reductions.py [full_path_to_raw_data] [full_path_to_reduced_data]

### Create wave maps for each of the arc frames
    python wht_measure_arc.py [arc file]

### Apply the wavemaps to object frames
    python wht_calibrate_objects.py [full_path_to_reduced_data]
    
### Rectify all the images in the current directory
    python apply_wavelength.py

