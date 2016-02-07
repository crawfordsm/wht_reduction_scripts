import sys
import numpy as np
from astropy.io import fits

def rectify(hdu):
    """Giving an object with a wave map, rectify the image
       and create a wcs for the image
    """
    data = hdu[0].data
    wmap = hdu[1].data

    i_c = int(data.shape[0]/2.0)
    warr = wmap[i_c,:]
    xarr = np.arange(len(warr))
    dw = (warr.max()-warr.min())/len(xarr)
    warr = np.arange(warr.min(), warr.max(), dw)
    for i in range(data.shape[0]):
        data[i,:] = np.interp(warr, wmap[i,:], data[i,:])
    hdu[0].data = data
    hdu[0].header['CTYPE1'] = 'LAMBDA'
    hdu[0].header['CTYPE2'] = 'PIXEL'
    hdu[0].header['CD1_1'] = dw
    hdu[0].header['CD2_1'] = 0.0
    hdu[0].header['CD1_2'] = 0.0
    hdu[0].header['CD2_2'] = 1.0 
    hdu[0].header['CRPIX1'] = 0.0
    hdu[0].header['CRPIX2'] = 0.0
    hdu[0].header['CRVAL1'] = warr.min()
    hdu[0].header['CRVAL2'] = 0.0
    hdu[0].header['CDELT1'] = 1.0
    hdu[0].header['CDELT2'] = 1.0
    hdu[0].header['DC-FLAG'] = 0

    return hdu


if __name__=='__main__':
   #filename=sys.argv[1]
   import glob
   import ccdproc
   for filename in glob.glob('w_r*fit'):
      hdu = fits.open(filename)
      hdu[0].data, crmask = ccdproc.cosmicray_lacosmic(hdu[0].data)
      hdu = rectify(hdu)
      hdu.writeto(filename.replace('w_', 't_'), clobber=True)
      

