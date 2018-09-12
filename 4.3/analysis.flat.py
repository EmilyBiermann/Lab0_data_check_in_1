"""

Emily Biermann
AST 443
Lab0
Data Analysis
4.3 Imaging Flat-Fields

"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

pwd = Path("/home/ebiermann/AST443/Lab0/Lab0_data")

# Flat Fields
subfolder = "3.2_CCDFlats"
filenames = ['3.2_flat_real.00000000.FIT', '3.2_flat_real.00000001.FIT',
             '3.2_flat_real.00000002.FIT', '3.2_flat_real.00000003.FIT',
             '3.2_flat_real.00000004.FIT', '3.2_flat_real.00000005.FIT',
             '3.2_flat_real.00000006.FIT', '3.2_flat_real.00000007.FIT',
             '3.2_flat_real.00000008.FIT', '3.2_flat_real.00000009.FIT' ]

### 1 ###

# get master bias from 4.1
pwd2 = Path("/home/ebiermann/AST443/Lab0/Scripts/4.1")
masterbias = np.load(pwd2 / 'm10_masterbias.npy')

# open files and extract images
allimages = []
N=0 # number of frames
for file in filenames:
    N=N+1
    fname = pwd / subfolder / file
    hdu = fits.open(fname)
    allimages.append(hdu[0].data)

# Create masterflat from mean combine
masterflat = np.mean(allimages) - masterbias

### 2 ###

# save masterflat for use in ds9
#masterflat_FIT = fits.PrimaryHDU(masterflat)
#masterflat_FIT.writeto('masterflat.fits')

### 3 ###

# Flatten masterflat
countval = masterflat.flatten()

# make histogram of masterflat to identify dead pixels
plt.figure(1)
plt.title(r'Master Flat')
plt.xlabel(r'Number of Counts')
plt.ylabel(r'Numberof Pixels')
plt.hist(countval,bins=100)
plt.yscale('log')
plt.savefig('masterflat.png',format='png',dpi=1000,bbox_inches='tight')

### 4 ###

# load data
data = np.loadtxt('sensitivity_data.txt',skiprows=1)
x = data[0,:] - data[0,0]
y = data[1,:] - data[0,1]
rel_brightness = data[2,:] - data[2,0]

# calculate distance from center
d = np.sqrt(x**2+y**2)

# plot sensitivity vs distance
plt.figure(2)
plt.title(r'Camera Sensitivity')
plt.xlabel(r'Distance from Center')
plt.savefig('sensitivity_plot.png',format='png',dpi=1000,bbox_inches='tight')


plt.errorbar(d,rel_brightness,fmt='x')

### 5 ###

# open 90deg rotation image
fname = pwd / subfolder / '3.2_flat_real_rotated.00000010.FIT'
hdu = fits.open(fname)
rotflat = hdu[0].data

# subtract bias
rotflat = rotflat - masterbias

# save rotated flat image
#masterflat_FIT = fits.PrimaryHDU(rotflat)
#masterflat_FIT.writeto('rotflat.fits')


