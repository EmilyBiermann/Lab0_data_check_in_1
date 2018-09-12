"""

Emily Biermann
AST 443
Lab0
Data Analysis
4.4 Bad Pixel Map

Makes a binary bad pixel map for a FIT file

"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

# Path and file name
pwd = Path("/home/ebiermann/AST443/Lab0/Lab0_data/3.2_CCDFlats")
fname = pwd / '3.2_flat_real.00000000.FIT'

# Open image
hdu = fits.open(fname)

# Get biased
# get master bias from 4.1
pwd2 = Path("/home/ebiermann/AST443/Lab0/Scripts/4.1")
masterbias = np.load(pwd2 / 'm10_masterbias.npy')

# extract data and subtract bias
imagedata = hdu[0].data - masterbias

# define cuts

imagedata_flat = imagedata.flatten()
imagedata_stats=stats.sigmaclip(imagedata_flat,5.0,5.0)

# use beyond 5sigma as hotcut
hotcut = imagedata_stats[2]

# use below 5sigma as dead cut
deadcut = imagedata_stats[1]

bad_pixel_map = np.zeros_like(imagedata)
count=0
for i in range(len(imagedata[0])):
    for j in range(len(imagedata[1])):
        if imagedata[i,j] < deadcut or imagedata[i,j] > hotcut:
            bad_pixel_map[i,j] = 1
            count=count+1
        else:
            continue

# save bad_pixel_map as fit file
badpixelmap_FIT = fits.PrimaryHDU(bad_pixel_map)
badpixelmap_FIT.writeto('bad_pixel_map.fits')
                
