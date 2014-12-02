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
parser.add_argument('-n','--normalise', help='normalises the profiles to the peak', action='store_true')
parser.add_argument('-l','--listbadmjds', help='file with list of observation dates which are to be excluded from the dataset')

args = parser.parse_args()
print 'Read arguments',
filename = args.filename
pulsar = args.pulsar
print 'PSR',pulsar
outlierthreshold = args.removeoutliers
binstartzoom = []
binendzoom = []

rawdata = np.loadtxt(filename)


# Remove any observation on the mjds listed in a txt file

if args.listbadmjds:

    bad_mjd_file = args.listbadmjds
    bad_mjds_index = []
    bad_mjds = np.loadtxt(bad_mjd_file)
    bad_mjds = np.atleast_1d(bad_mjds)
    for b in bad_mjds:
        for i in range(rawdata.shape[1]):
            if abs(rawdata[rawdata.shape[0]-1,i]-b) < 1.0:
                bad_mjds_index.append(i)

    rawdata = np.delete(rawdata,bad_mjds_index,1)



mjd = rawdata[rawdata.shape[0]-1,:]
Tobs = rawdata[rawdata.shape[0]-2,:]
data = rawdata[0:rawdata.shape[0]-2,:]
bins = data.shape[0]

if not (os.path.exists('./{0}/'.format(pulsar))):
    os.mkdir('./{0}/'.format(pulsar))  

if (args.normalise):
    print "Pulse profiles WILL be normalised"
else:
    print "Pulse profiles will NOT be normalised"

# Make plots of originals if needed
if (args.originalplots):
    dir='raw_profiles'
    Vf.makeplots(pulsar,data,mjd,dir,bins)

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

suffix = ''

# Normalise data to peak if flag -n used. Peak should already be aligned and at nbins/4:

if (args.normalise):
    aligned_data = Vf.norm_to_peak(aligned_data,aligned_data.shape[0]/4-1)

    for i in range(aligned_data.shape[0]):
        originaltemplate[i]=np.median(aligned_data[i,:])

    suffix = '_norm'

# Find pulse regions on template
rmstemplate = np.median(rmsperepoch) / np.sqrt(baselineremoved.shape[1])
peaks = 1
regioncounter = 0
peakoriginal = np.max(originaltemplate)

while peaks != 0 and regioncounter < 5:
    bs, be, peaks, cuttemplate = Vf.binstartend(template, peakoriginal, rmstemplate)
    binstartzoom.append(bs)
    binendzoom.append(be)
    template = cuttemplate
    regioncounter += 1
# last attempt has failed, hence we are here, so:
regioncounter -= 1

if pulsar == 'J0738-4042':
    binstartzoom = [138]
    binendzoom = [357]
if pulsar == 'J0742-2822':
    binstartzoom = [222]
    binendzoom = [307]
if pulsar == 'J0842-4851':
    binstartzoom = [212, 721]
    binendzoom = [297,806]
if pulsar == 'J0908-4913':
    binstartzoom = [202,708]
    binendzoom = [292,798]
if pulsar == 'J0940-5428':
    binstartzoom = [16]
    binendzoom = [46]
if pulsar == 'J1105-6107':
    binstartzoom = [181]
    binendzoom = [292]
if pulsar == 'J1302-6350':
    binstartzoom = [6, 64]
    binendzoom = [49, 107]
if pulsar == 'J1359-6038':
    binstartzoom = [202]
    binendzoom = [307]
if pulsar == 'J1600-5044':
    binstartzoom = [217]
    binendzoom = [318]
if pulsar == 'J1602-5100':
    binstartzoom = [227]
    binendzoom = [292]
if pulsar == 'J1730-3350':
    binstartzoom = [18]
    binendzoom = [61]
if pulsar == 'J1740-3015':
    binstartzoom = [227]
    binendzoom = [287]
if pulsar == 'J1757-2421':
    binstartzoom = [165]
    binendzoom = [306]
if pulsar == 'J1825-0935':
    binstartzoom = [180,760]
    binendzoom = [290,810]
if pulsar == 'J1830-1059':
    binstartzoom = [227]
    binendzoom = [282]
if pulsar == 'J1845-0743':
    binstartzoom = [151]
    binendzoom = [327]

print 'Found ', regioncounter, ' regions'
left = np.array(binstartzoom)
right = np.array(binendzoom)
binline = np.zeros(3)
for i in range(regioncounter):
    outputarray = aligned_data[left[i]:right[i],:]
    outputfile = '{0}/zoomed_{0}_bins_{1}-{2}{3}.txt' .format(pulsar,left[i],right[i],suffix)
    np.savetxt(outputfile, outputarray)
    if (args.diagnosticplots):
        plt.imshow(outputarray,aspect = 'auto')
        plt.colorbar(orientation="horizontal") 
        plt.savefig('{0}/zoomed_{0}_bins_{1}-{2}{3}.png' .format(pulsar,left[i],right[i],suffix))
        plt.clf()
    # Write info file in created directory
    binline[0] = left[i]
    binline[1] = right[i]
    binline[2] = aligned_data.shape[0]
    outputfile = '{0}/zoomed_{0}_bins_{1}-{2}{3}.dat' .format(pulsar,left[i],right[i],suffix)
    np.savetxt(outputfile, binline)

    # If profiles deviate far from the median profile, they are flagged and saved in a folder called 'flagged_profiles'
    if (args.normalise):

        resamp_for_flagged = scisig.resample(aligned_data,aligned_data.shape[0]/8)
        resamp_template = scisig.resample(originaltemplate,aligned_data.shape[0]/8)
        if i == 0:
            sum_diff = np.zeros((aligned_data.shape[1]))

            for k in range(aligned_data.shape[1]):
                for j in range(resamp_for_flagged.shape[0]):
                    sum_diff[k]+=abs(resamp_for_flagged[j,k]-resamp_template[j])

            os.mkdir('./{0}/flagged_profiles'.format(pulsar))

            flagged = 0

        print "SHAPE OF DATA AND MJD",aligned_data.shape[1],mjdout.shape[0]

        for k in range(aligned_data.shape[1]):
            if sum_diff[k] > np.mean(sum_diff[:])+1.5*np.std(sum_diff[:]):
                plt.plot(aligned_data[left[i]:right[i],k],'b-')
                plt.plot(originaltemplate[left[i]:right[i]],'r-')
                plt.savefig('./{0}/flagged_profiles/{1}_zoom_{2}.png' .format(pulsar,mjdout[k],i))
                plt.clf()
#                plt.plot(resamp_for_flagged[left[i]/8:right[i]/8,k],'b-')
#                plt.plot(resamp_template[left[i]/8:right[i]/8],'r-')
#                plt.savefig('./{0}/flagged_profiles/{1}_zoom_binned_{2}.png' .format(pulsar,mjdout[k],i))

#                plt.clf()
                if i == 0:
                    flagged+=1
                    plt.plot(aligned_data[:,k],'b-')
                    plt.plot(originaltemplate,'r-')
                    plt.savefig('./{0}/flagged_profiles/{1}.png' .format(pulsar,mjdout[k]))
                    plt.clf()
#                    plt.plot(resamp_for_flagged[:,k],'b-')
#                    plt.plot(resamp_template,'r-')
#                    plt.savefig('./{0}/flagged_profiles/{1}_binned.png' .format(pulsar,mjdout[k]))
#                    plt.clf()

        if i == 0:
            print "Percentage of profiles flagged as unusual:",np.round(float(flagged)/aligned_data.shape[1]*100.0,2),"%",flagged,aligned_data.shape[1]

if (args.goodprofiles):
    dir='good_profiles{0}' .format(suffix)

    if regioncounter > 1:
        # Make interpulse plots of good profiles if needed
        Vf.goodplots_ip(pulsar,aligned_data[left[0]:right[0],:], aligned_data[left[1]:right[1],:],mjdout,dir,bins,left[1],originaltemplate[left[0]:right[0]], originaltemplate[left[1]:right[1]],-0.5*np.std(aligned_data[left[0]:right[0],:]),np.max(aligned_data[left[0]:right[0],:]),peakindex=bins/4-left[0])

    else:
        # Make plots of good profiles if needed
        Vf.makeplots(pulsar,aligned_data[left[0]:right[0],:],mjdout,dir,bins,template=originaltemplate[left[0]:right[0]],yllim=-0.1*np.max(originaltemplate),yulim=np.max(aligned_data[left[0]:right[0],:]),peakindex=bins/4-left[0])

#to make plots for the paper
#        Vf.forpaper(pulsar,aligned_data[left[0]:right[0],162],aligned_data[left[0]:right[0],40],originaltemplate[left[0]:right[0]],yllim=-0.1*np.max(originaltemplate),yulim=np.max(aligned_data[left[0]:right[0],:]))

# Make plots of removed profiles if needed
if (args.badprofiles):
    dir='removed_profiles{0}' .format(suffix)
    Vf.makeplots(pulsar,removedprofiles,mjdremoved,dir,bins,template=originaltemplate,yllim=-0.1*np.max(originaltemplate),yulim=np.max(aligned_data[left[0]:right[0],:]))

outputfile = '{0}/mjd{1}.txt' .format(pulsar,suffix)
np.savetxt(outputfile, mjdout)
outputfile = '{0}/mjdremoved{1}.txt' .format(pulsar,suffix)
np.savetxt(outputfile, mjdremoved)

if (args.diagnosticplots):
    plt.imshow(data,aspect = 'auto')
    plt.colorbar(orientation="horizontal")
    plt.savefig('./{0}/{0}_rawdata.png' .format(pulsar))
    plt.imshow(baselineremoved,aspect = 'auto')
    plt.savefig('./{0}/{0}_baselined.png' .format(pulsar))
