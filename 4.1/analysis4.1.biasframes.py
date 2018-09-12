#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 12:39:53 2018

Emily Biermann
AST 443
Lab0
Data Analysis
4.1 Bias Frames

"""

# 1 #

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

pwd = Path("/home/ebiermann/AST443/Lab0/Lab0_data")

"""

# -10degC Bias Frames
subfolder = "3.1_m10_BIAS"
filenames = [ '3.1_m10.00000000.BIAS.FIT', '3.1_m10.00000005.BIAS.FIT',
              '3.1_m10.00000001.BIAS.FIT', '3.1_m10.00000006.BIAS.FIT',
              '3.1_m10.00000002.BIAS.FIT', '3.1_m10.00000007.BIAS.FIT',
              '3.1_m10.00000003.BIAS.FIT', '3.1_m10.00000008.BIAS.FIT',
              '3.1_m10.00000004.BIAS.FIT', '3.1_m10.00000009.BIAS.FIT'  ]
fig1 = '4.1_m10_BIAS_hist.png'
fig2 = '4.1_m10_BIAS_histfit.png'
mastername = 'm10_masterbias'

"""

# +10degC Bias Frames
subfolder = "3.1_p10_BIAS"
filenames = [ '3.1_p10.00000000.BIAS.FIT', '3.1_p10.00000005.BIAS.FIT',
              '3.1_p10.00000001.BIAS.FIT', '3.1_p10.00000006.BIAS.FIT',
              '3.1_p10.00000002.BIAS.FIT', '3.1_p10.00000007.BIAS.FIT',
              '3.1_p10.00000003.BIAS.FIT', '3.1_p10.00000008.BIAS.FIT',
              '3.1_p10.00000004.BIAS.FIT', '3.1_p10.00000009.BIAS.FIT'  ]
fig1 = '4.1_p10_BIAS_hist.png'
fig2 = '4.1_p10_BIAS_histfit.png'
mastername = 'p10_masterbias'

#"""

# open one Bias Frame
fname = pwd / subfolder / filenames[0]
bias1 = fits.open(fname)

# extract and flatten image data
imagedata = bias1[0].data
countval = imagedata.flatten()

# plot histogram of count values
plt.figure(1)
plt.title(r'Number of Counts per Pixel')
plt.xlabel(r'Number of Counts')
plt.ylabel(r'Numberof Pixels')
plt.hist(countval, range=[900,1300], bins=100)
plt.yscale('log')
plt.savefig(fig1,format='png',dpi=1000,bbox_inches='tight')

# Calculate and print Statistics 
mean1 = np.mean(countval)
print('mean =',mean1)

print('mode =',stats.mode(countval)[0][0])

print('median =',np.median(countval))

std1 = np.std(countval)
print('std =',std1)

# Calculate histogram parameters
cmin=900
cmax=1200
nbins=100
nrm=(cmax-cmin)/nbins*len(countval[(countval>=cmin) & (countval<=cmax)])
clipmin=cmin
clipmax=1100
clipval=countval[(countval>=clipmin) & (countval<=clipmax)]
percent_clip=1-len(clipval)/len(countval)
print('fraction of rejected pixels =',percent_clip)

mu=np.mean(clipval)
sig=np.std(clipval)
mode=stats.mode(clipval)[0][0]

xarray=np.linspace(cmin,cmax,nbins*10)
yarray=nrm*norm.pdf(xarray,loc=mu,scale=sig)

# plot histogram with fitted gaussian
plt.figure(2)
plt.title(r'Number of Counts per Pixel with Fitted Histogram')
plt.xlabel(r'Number of Counts')
plt.ylabel(r'Numberof Pixels')
plt.hist(countval,range=[cmin,cmax],bins=nbins)
plt.yscale('log')
plt.ylim([0.1,1e6])
plt.plot(xarray,yarray,color='red',linewidth=2.0)
plt.axvline(x=mode,color='green',linewidth=2.0)
plt.savefig(fig2,format='png',dpi=1000,bbox_inches='tight')

# 2 #

# Find Gain
header = bias1[0].header
gain = header['EGAIN'] # N_elec/N_counts
print('gain =',gain)

# Calculate read noise
readnoise = sig*gain # convert to units of electrons
print('read noise =', readnoise, 'e-rms')
manread = 15.75
print('man read noise =', manread, 'e-rms')
    # man: 15e-rms
print('percent diff = ', ((readnoise-manread)/manread)*100.0)

### 3 ###

# allocate array for images
allimages = []

# open files and extract images
N=0 # number of frames
for file in filenames:
    N=N+1
    fname = pwd / subfolder / file
    hdu = fits.open(fname)
    allimages.append(hdu[0].data)

# create master bias from 10 images
masterbias = np.mean(allimages,axis=0)

allcounts = masterbias.flatten()

# calculate and print mean, std
mean_all = np.mean(allcounts,axis=0)
std_all = np.std(allcounts,axis=0)
print('Total Bias mean =', mean_all)
print('Total Bias std =', std_all)
    
# std dec. factor
std_decfac = (std1-std_all)/std1
print('std decr. fac =', std_decfac)

# expected dec. factor
print('1/sqrt(N) =', 1/np.sqrt(N))

np.save(mastername,masterbias)




