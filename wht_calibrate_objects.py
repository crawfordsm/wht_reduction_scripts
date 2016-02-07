import sys
import os
import numpy as np
  
from astropy.io import fits
from astropy import units as u

import ccdproc 
from ccdproc import CCDData
from datetime import datetime

from ccdproc import ImageFileCollection

if len(sys.argv)!=2: 
   print('Usage:\npython wht_calibrate_objects.py [full_path_to_reduced_data]\n')
   exit()

indir = sys.argv[1]

os.chdir(indir)

#change this to point to your raw data directory
ic1 = ImageFileCollection(indir.replace('sci2', '20160115'))

#create an array of dates
arm = 'Red arm'
file_list = []
date_list=[]
for hdu, fname in ic1.hdus(obstype='Arc', isiarm=arm, return_fname=True):
    if os.path.isfile('w_arc_'+os.path.basename(fname)):
        d = hdu.header['DATE-OBS'] + ' ' + hdu.header['UT']
        d = datetime.strptime(d, '%Y-%m-%d %H:%M:%S.%f')
        file_list.append(fname)
        date_list.append(d)
date_arr = np.array(date_list)

#reduce the object frames
for filename in ic1.files_filtered(obstype='TARGET', isiarm=arm):
    hdu = fits.open('obj_' + filename)

    #find the closest two arcs in time
    #open them, based on the time difference 
    #average their frames, and then use that as the 
    #wave map for the object file

    date = hdu[0].header['DATE-OBS'] + ' ' + hdu[0].header['UT']
    date = datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f')
    d = abs(date_arr - date)
    i = d.argmin()
    if date_arr[i] < date:
       j = i+1
    else:
       j = i-1
    w1 = 1 - 1.0*d[i].seconds/(d[i].seconds+d[j].seconds)
    w2 = 1 - 1.0*d[j].seconds/(d[i].seconds+d[j].seconds)
    wdata1 = fits.open('w_arc_' + os.path.basename(file_list[i]))['WAV'].data
    wdata2 = fits.open('w_arc_' + os.path.basename(file_list[j]))['WAV'].data
    wdata = w1*wdata1 + w2*wdata2

    hdu.append(fits.ImageHDU(data=wdata, name='WAV'))
    hdu.writeto('w_' + os.path.basename(filename), clobber=True)

