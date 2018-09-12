"""

Emily Biermann
AST 443
Lab0
Data Analysis
4.2 Dark Frames (-10degC)

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

# -10degC Dark Frames
subfolder = "3.1_m10_DARK_1minexpo"
filenames = ['3.1_m10.00000017.DARK.FIT',  '3.1_m10.00000022.DARK.FIT',
             '3.1_m10.00000018.DARK.FIT',  '3.1_m10.00000023.DARK.FIT',
             '3.1_m10.00000019.DARK.FIT',  '3.1_m10.00000024.DARK.FIT',
             '3.1_m10.00000020.DARK.FIT',  '3.1_m10.00000025.DARK.FIT',
             '3.1_m10.00000021.DARK.FIT',  '3.1_m10.00000026.DARK.FIT' ]

# -10degC Darks with Varying Exposure
subfolder2 = "3.1_m10_DARK_varyexpo"
filenames2 =  ['3.1_m10.00000010.DARK.FIT',  '3.1_m10.00000011.DARK.FIT',
               '3.1_m10.00000012.DARK.FIT',  '3.1_m10.00000013.DARK.FIT',
               '3.1_m10.00000014.DARK.FIT',  '3.1_m10.00000015.DARK.FIT',
               '3.1_m10.00000016.DARK.FIT' ]

fig1 = '4.2.1_m10_DARK_hist.png'
title1 = r'Dark Frames with 1min Exopsure at $-10^{\circ}C$'

fig2 = '4.2.1_m10_darkcurrenttest.png'
title2 = r'Dark Current Image Test'

fig3 = '4.2.1_m10_countsvsexpo.png'
title3 = r'Typical Counts vs Exposure Time ($-10^{\circ}C$)'

masterbias_name = 'm10_masterbias.npy'

# open files and extract images
allimages = []
N=0 # number of frames
for file in filenames:
    N=N+1
    fname = pwd / subfolder / file
    hdu = fits.open(fname)
    allimages.append(hdu[0].data)

### 1 ###

# median combine allimages to create master dark frame
masterdark = np.median(allimages,axis=0)
countval = masterdark.flatten()

# Make histogram of counts
plt.figure(1)
plt.title(title1)
plt.xlabel(r'Number of Counts')
plt.ylabel(r'Numberof Pixels')
plt.hist(countval, range=[950,5000], bins=100)
plt.yscale('log')
plt.savefig(fig1,format='png',dpi=1000,bbox_inches='tight')

# Calculate and print Statistics 
mean = np.mean(countval)
print('mean =',mean)
mode = stats.mode(countval)[0][0]
print('mode =',mode)
median = np.median(countval)
print('median =',median)
std = np.std(countval)
print('std =',std)

# Define Cut
cmin=950
cmax=4000 
clipval_stats=stats.sigmaclip(countval,5.0,5.0)
clipval=clipval_stats[0]
frac_clip=1-len(clipval)/len(countval)
print('Fraction of rejected pixels =', frac_clip)

# Calculate and print Statistics from clipped set
cmean = np.mean(clipval)
print('clip mean =',cmean)
cmode = stats.mode(clipval)[0][0]
print('clip mode =',cmode)
cmedian = np.median(clipval)
print('clip median =',cmedian)
cstd = np.std(clipval)
print('clip std =',cstd)

# Hot pixels
hotcut = clipval_stats[2]
hotpix = countval[countval>hotcut]
frac_hot = len(hotpix)/len(countval)
print('Hot Cut =', hotcut)
print('Fraction of Hot Pixels =', frac_hot)

### 2 ###

# get master bias from 4.1
pwd2 = Path("/home/ebiermann/AST443/Lab0/Scripts/4.1")
masterbias = np.load(pwd2 / masterbias_name)

# initialize array to hold all darkcurrent images
darkcurrent_all = []

# compute dark current by subtracting masterbias from dark frames
for dark in allimages:
    darkcurrent = dark - masterbias
    darkcurrent = np.asarray(darkcurrent)
    darkcurrent_val = darkcurrent.flatten()
    
# make histogram of one to determine method ie mean, median, mode
darkcurrent_val = darkcurrent.flatten()
plt.figure(2)
plt.title(title2)
plt.hist(darkcurrent_val,bins=100)
plt.yscale('log')
plt.savefig(fig2,format='png',dpi=1000,bbox_inches='tight')
    # Let's use mode, too many outliers for mean or median

# Found mode of all darkcurrent images


### 3 ###

# open files and extract images, exposure time
allimages = []
exptime = []
N=0 # number of frames
for file in filenames2:
    N=N+1
    fname = pwd / subfolder2 / file
    hdu = fits.open(fname)
    exptime.append(hdu[0].header['EXPTIME'])
    allimages.append(hdu[0].data)

# initialize array to hold mode of each image, uncertainty, exptime
darkcurrent_mode = []
darkcurrent_sig = []

# compute dark current by subtracting masterbias from dark frames
for dark in allimages:
    darkcurrent = dark - masterbias
    darkcurrent = np.asarray(darkcurrent)
    darkcurrent_val = darkcurrent.flatten()
    darkcurrent_mode.append(stats.mode(darkcurrent_val)[0][0])
    # using sqrt(mean) bc poisson dist.
    darkcurrent_sig.append(np.sqrt(np.mean(darkcurrent_val)))

# perform linear regression
linearfit = np.polynomial.polynomial.polyfit(exptime,darkcurrent_mode, deg=1)
x=np.linspace(0,350,10**3)
def line(x,b,m):
    return b+m*x
lin_reg = line(x,*linearfit)

# plot data and fit
plt.figure(3)
plt.title(title3)
plt.errorbar(exptime,darkcurrent_mode,yerr=darkcurrent_sig,fmt='x',label='data')
plt.plot(x,lin_reg, label='linear fit')
plt.xlabel(r'Exposure Time [s]')
plt.ylabel(r'Mode Dark Current Counts [e-]')
plt.legend(loc='best')
plt.savefig(fig3,format='png',dpi=1000,bbox_inches='tight')

# print slope
print(' ')
print('darkcurrent/pixel/sec =', linearfit[1])
    
