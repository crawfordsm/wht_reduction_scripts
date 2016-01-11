import sys
import numpy as np
import argparse

from ccdproc import CCDData
from pyhrs import red_process

import specutils
from astropy import units as u

from astropy import modeling as mod
from astropy.io import fits

import pylab as pl

import specreduce
from specreduce.models import ISISModel

from specreduce.interidentify import InterIdentify
from specreduce import spectools as st
from specreduce import WavelengthSolution


#Xe_spec = specutils.io.read_ascii.read_ascii_spectrum1d('Xe.salt', dispersion_unit=u.angstrom)


parser = argparse.ArgumentParser(description='Reduce arc data from WHT ISIS')
parser.add_argument('arcfile', help='an arcfile to wavelength calibrate')
parser.add_argument('--list', dest='arc_ref', default='CuArNe.txt', 
                   help='File with list of known arc lines')

args = parser.parse_args()
arcfile=args.arcfile
print arcfile, type(arcfile)
arc = fits.open(arcfile)

#read in arc lines
arc_ref=args.arc_ref
print arc_ref
slines, sfluxes = np.loadtxt(arc_ref, usecols=(0,1), unpack=True)

function='poly'
order=3
rstep=1
nrows=1
mdiff=20
wdiff=3
thresh=3
niter=5
dc=3
ndstep=50
dsigma=5
method='Zeropoint'
res=6.0
dres=0.1
filename=None
smooth=0
inter=True
subback=0
textcolor='black'
log = None




#set up the model for the spectrograph
data = arc[0].data
header = arc[0].header
xarr = np.arange(data.shape[1])

instrume = header['INSTRUME'].strip()

#

grating = header['ISIGRAT'].strip()
arm = header['ISIARM']
slit = header['ISISLITW']
# i don't know exactly how the ISITHETA relats
# to the grating angle alpha so we are actually 
# just going to use the information in the header
import math
wc =  header['CENWAVE'] * u.angstrom
s =  1.0/header['LINESMM'] * u.mm
gratang = 180 / math.pi * math.asin((wc/s).cgs/2.0 * 1.0)
print gratang
if arm.count('Red'):
   grang = 23.1824 # header['ISITHETA']/1000.0) 
else:
   grang = 17.1 # header['ISITHETA']/1000.0) 
xbin = 1
ybin = 1
xpos = 0 #-0.2666
ypos = 0 #0.0117
print instrume, grating, header['ISITHETA']/1000.0, grang, arm, slit, header['CENWAVE'], header['DISPERSI']

rss = ISISModel.ISISModel(grating_name=grating.strip(), gratang=grang,  camera_name=arm,
                                            slit=slit, xbin=xbin, ybin=ybin,
                                            xpos=xpos, ypos=ypos)
data = arc[0].data
xarr = np.arange(data.shape[1])
warr = rss.get_wavelength(xarr) * u.mm
warr = warr.to(u.angstrom)
print int(0.5*len(xarr)),  warr[int(0.5*len(xarr))], rss.calc_centralwavelength()
lmax = data[500].max()
swarr, sfarr = st.makeartificial(slines, sfluxes, lmax, res, dres)
pl.plot(warr, data[500])
mask = (swarr > warr.value.min()) * ( swarr < warr.value.max())
#pl.plot(swarr[mask], sfarr[mask])
#pmel.show()
#nws = fit_ws(ws_init, xarr, warr)


ws_init = mod.models.Legendre1D(3)
ws_init.domain = [xarr.min(), xarr.max()]
ws = WavelengthSolution.WavelengthSolution(xarr, warr.value, ws_init)
ws.fit()

istart = data.shape[0]/2.0
smask = (slines > warr.value.min()-10) * (slines < warr.value.max()+10)
iws = InterIdentify(xarr, data, slines[smask], sfluxes[smask], ws, mdiff=mdiff, rstep=rstep,
              function=function, order=order, sigma=thresh, niter=niter, wdiff=wdiff,
              res=res, dres=dres, dc=dc, ndstep=ndstep, istart=istart,
              method=method, smooth=smooth, filename=filename,
              subback=subback, textcolor=textcolor, log=log, verbose=True)

import pickle
name = sys.argv[1].replace('fit', 'pkl')
pickle.dump(iws, open(name, 'w'))

ws_init = mod.models.Legendre1D(3)
ws_init.domain = [xarr.min(), xarr.max()]
ws = WavelengthSolution.WavelengthSolution(xarr, xarr, ws_init)
ws.fit()
istart = int(0.5 * len(data))
aws = st.arc_straighten(data, istart, ws, rstep=1)

data = st.wave_map(data, aws)
k = iws.keys()[0]
for i in range(data.shape[0]):
    data[i,:] = iws[k](data[i,:])

arc.append(fits.ImageHDU(data=wave_map, header=arc['SCI'].header, name='WAV'))
arc.writeto('w_' + os.path.basename(arcfile), clobber=True)

exit()


