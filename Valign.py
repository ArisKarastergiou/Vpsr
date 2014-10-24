#!/usr/bin/python
# Reads in 2d profile data(epoch, bin), aligns, normalizes, rejects outliers and writes out

import argparse
import pylab as pb
pb.ion()
import numpy as np
import GPy
import math
import matplotlib
import matplotlib.pyplot as plt
import scipy as sc
import scipy.signal as scisig
import sys
import Vfunctions as Vf
import os 

# Read command line arguments
parser = argparse.ArgumentParser(description='Pulsar profile alignment for variability studies')
parser.add_argument('-f','--filename', help='File 2d data set', required=True)
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-o','--originalplots', help='make plots of raw data', action='store_true')
parser.add_argument('-b','--badprofiles', help='make plots of bad, removed profiles', action='store_true')
parser.add_argument('-g','--goodprofiles', help='make plots of good profiles', action='store_true')
parser.add_argument('-d','--diagnosticplots', help='make image plots', action='store_true')
parser.add_argument('-a','--autozoom', help='number of regions to autozoom', action='store_true')
parser.add_argument('-r','--removeoutliers', help='factor of deviation from median baseline rms to remove outliers ', type=float, required='True')

args = parser.parse_args()

print 'Read arguments'
filename = args.filename
pulsar = args.pulsar
outlierthreshold = args.removeoutliers
binstartzoom = []
binendzoom = []

rawdata = np.loadtxt(filename)


mjd = rawdata[rawdata.shape[0]-1,:]
Tobs = rawdata[rawdata.shape[0]-2,:]
data = rawdata[0:rawdata.shape[0]-2,:]
bins = data.shape[0]

if not (os.path.exists('./{0}/'.format(pulsar))):
    os.mkdir('./{0}/'.format(pulsar))  

# Make plots of originals if needed
if (args.originalplots):
    dir='raw_profiles'
    Vf.makeplots(pulsar,data,mjd,dir)

# Remove baseline and outliers if requested
baselineremoved, removedprofiles, rmsperepoch, outlierlist, inlierlist = Vf.removebaseline(data, outlierthreshold)
mjdout = np.delete(mjd,outlierlist)
mjdremoved = np.delete(mjd,inlierlist)

print 'Baseline and outliers removed'
print 'Remaining array shape:', baselineremoved.shape
brightestprofile = Vf.findbrightestprofile(baselineremoved,rmsperepoch)
print 'Brightest profile: ', brightestprofile

# resample if brightest profile is less than 20 sigma peak
brightprofpeak = np.max(baselineremoved[:,brightestprofile])
brightprofrms = rmsperepoch[brightestprofile]
if brightprofpeak/brightprofrms < 20 :
    resampled = scisig.resample(baselineremoved,bins/8)
    print 'resampling 8'
else:
    resampled = baselineremoved
baselineremoved = resampled


# Align data by the following way: x-correlate profiles with brightest, shift then x-correlate with average

aligned_data, template = Vf.aligndata(baselineremoved, brightestprofile, pulsar)

originaltemplate = np.copy(template)

# Find pulse regions on template
rmstemplate = np.mean(rmsperepoch) / np.sqrt(baselineremoved.shape[1])
peaks = 1
regioncounter = 0

while peaks != 0:
    bs, be, peaks, cuttemplate = Vf.binstartend(template, rmstemplate)
    binstartzoom.append(bs)
    binendzoom.append(be)
    template = cuttemplate
    regioncounter += 1
# last attempt has failed, hence we are here, so:
regioncounter -= 1

print 'Found ', regioncounter, ' regions'
left = np.array(binstartzoom)
right = np.array(binendzoom)
binline = np.zeros(3)
for i in range(regioncounter):
    outputarray = aligned_data[left[i]:right[i],:]
    outputfile = '{0}/zoomed_{0}_bins_{1}-{2}.txt' .format(pulsar,left[i],right[i])
    np.savetxt(outputfile, outputarray)
    if (args.diagnosticplots):
        plt.imshow(outputarray,aspect = 'auto')
        plt.colorbar(orientation="horizontal") 
        plt.savefig('{0}/zoomed_{0}_bins_{1}-{2}.png' .format(pulsar,left[i],right[i]))
        plt.clf()
    # Write info file in created directory
    binline[0] = left[i]
    binline[1] = right[i]
    binline[2] = aligned_data.shape[0]
    outputfile = '{0}/zoomed_{0}_bins_{1}-{2}.dat' .format(pulsar,left[i],right[i])
    np.savetxt(outputfile, binline)

# Make plots of good profiles if needed
if (args.goodprofiles):
    dir='good_profiles'
    Vf.makeplots(pulsar,aligned_data[binline[0]:binline[1],:],mjdout,dir,template=originaltemplate[binline[0]:binline[1]],yllim=-100,yulim=np.max(aligned_data[binline[0]:binline[1],:]),peakindex=bins/4-left[0])

# Make plots of removed profiles if needed
if (args.badprofiles):
    dir='removed_profiles'
    Vf.makeplots(pulsar,removedprofiles,mjdremoved,dir,template=originaltemplate,yllim=-100,yulim=np.max(aligned_data[binline[0]:binline[1],:]))

outputfile = '{0}/mjd.txt' .format(pulsar)
np.savetxt(outputfile, mjdout)

outputfile = '{0}/mjdremoved.txt' .format(pulsar)
np.savetxt(outputfile, mjdremoved)


if (args.diagnosticplots):
    plt.imshow(data,aspect = 'auto')
    plt.colorbar(orientation="horizontal")
    plt.savefig('./{0}/{0}_rawdata.png' .format(pulsar))
    plt.imshow(baselineremoved,aspect = 'auto')
    plt.savefig('./{0}/{0}_baselined.png' .format(pulsar))
