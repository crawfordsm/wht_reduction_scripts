import sys
import os
import numpy as np
  
from astropy.io import fits
from astropy import units as u

import ccdproc 
from ccdproc import CCDData

from ccdproc import ImageFileCollection

if len(sys.argv)!=3: 
   print('Usage:\npython wht_basic_rection.py [full_path_to_raw_data] [full_path_to_reduced_data]\n')
   exit()

indir = sys.argv[1]
outdir = sys.argv[2]

if not os.path.isdir(outdir): os.mkdir(outdir)
os.chdir(outdir)

#change this to point to your raw data directory
ic1 = ImageFileCollection(indir)

#create the bias frames
blue_bias_list = []
for filename in ic1.files_filtered(obstype='Bias', isiarm='Blue arm'):
    print ic1.location + filename
    ccd = CCDData.read(ic1.location + filename, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    blue_bias_list.append(ccd)
master_bias_blue = ccdproc.combine(blue_bias_list, method='median')
master_bias_blue.write('master_bias_blue.fits', clobber=True)

red_bias_list = []
for filename in ic1.files_filtered(obstype='Bias', isiarm='Red arm'):
    print ic1.location + filename
    ccd = CCDData.read(ic1.location + filename, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    red_bias_list.append(ccd)
master_bias_red = ccdproc.combine(red_bias_list, method='median')
master_bias_red.write('master_bias_red.fits', clobber=True)

#create the flat fields
red_flat_list = []
for filename in ic1.files_filtered(obstype='Flat', isiarm='Red arm'):
    ccd = CCDData.read(ic1.location + filename, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_red)
    red_flat_list.append(ccd)
master_flat_red = ccdproc.combine(red_flat_list, method='median')
master_flat_red.write('master_flat_red.fits', clobber=True)

blue_flat_list = []
for filename in ic1.files_filtered(obstype='Flat', isiarm='Blue arm'):
    ccd = CCDData.read(ic1.location + filename, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_blue)
    blue_flat_list.append(ccd)
master_flat_blue = ccdproc.combine(blue_flat_list, method='median')
master_flat_blue.write('master_flat_blue.fits', clobber=True)

#reduce the arc frames
for filename in ic1.files_filtered(obstype='Arc', isiarm='Blue arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_blue)
    ccd = ccdproc.flat_correct(ccd, master_flat_blue)
    ccd.data = ccd.data.T
    ccd.write('arc_'+filename, clobber=True)

red_flat_list = []
for filename in ic1.files_filtered(obstype='Arc', isiarm='Red arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_red)
    ccd = ccdproc.flat_correct(ccd, master_flat_red)
    ccd.data = ccd.data.T
    ccd.write('arc_'+filename, clobber=True)

#reduce the sky frames
for filename in ic1.files_filtered(obstype='Sky', isiarm='Blue arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_blue)
    ccd = ccdproc.flat_correct(ccd, master_flat_blue)
    ccd.data = ccd.data.T
    ccd.write('sky_'+filename, clobber=True)

for filename in ic1.files_filtered(obstype='Sky', isiarm='Red arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_red)
    ccd = ccdproc.flat_correct(ccd, master_flat_red)
    ccd.data = ccd.data.T
    ccd.write('sky_'+filename, clobber=True)

#reduce the object frames
for filename in ic1.files_filtered(obstype='TARGET', isiarm='Blue arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_blue)
    ccd = ccdproc.flat_correct(ccd, master_flat_blue)
    ccd.data = ccd.data.T
    ccd.write('obj_'+filename, clobber=True)

for filename in ic1.files_filtered(obstype='Target', isiarm='Red arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_red)
    ccd = ccdproc.flat_correct(ccd, master_flat_red)
    ccd.data = ccd.data.T
    ccd.write('obj_'+filename, clobber=True)
