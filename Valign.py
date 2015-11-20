#!/usr/bin/python
# Reads in 2d profile data(epoch, bin), aligns, normalizes, rejects outliers and writes out

import argparse
import pylab as pb
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy as sc
import scipy.signal as scisig
import sys
import Vfunctions as Vf
import os 
pb.ioff()

# Read command line arguments
parser = argparse.ArgumentParser(description='Pulsar profile alignment for variability studies')
parser.add_argument('-f','--filename', help='File 2d data set', required=True)
parser.add_argument('-p','--pulsar', help='Pulsar name', required=True)
parser.add_argument('-o','--originalplots', help='make plots of raw data', action='store_true')
parser.add_argument('-b','--badprofiles', help='make plots of bad, removed profiles', action='store_true')
parser.add_argument('-g','--goodprofiles', help='make plots of good profiles', action='store_true')
parser.add_argument('-d','--diagnosticplots', help='make image plots', action='store_true')
parser.add_argument('-r','--removeoutliers', help='factor of deviation from median baseline rms to remove outliers ', type=float, required='True')
parser.add_argument('-n','--normalise', help='normalises the profiles to the peak', action='store_true')
parser.add_argument('-l','--listbadmjds', help='file with list of observation dates which are to be excluded from the dataset')

args = parser.parse_args()
print 'Read arguments for',
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
    Vf.makeplots(pulsar,data,mjd,dir,bins,cal=1)

# Remove baseline and outliers if requested
baselineremoved, removedprofiles, rmsperepoch, outlierlist, inlierlist = Vf.removebaseline(data, outlierthreshold)
mjdout = np.delete(mjd,outlierlist)
mjdremoved = np.delete(mjd,inlierlist)

print 'Baseline and outliers removed'
print 'Removed MJD(s):',mjdremoved
print 'Remaining array shape (profiles, bins):', baselineremoved.shape
brightestprofile = Vf.findbrightestprofile(baselineremoved,rmsperepoch)

# resample if brightest profile is less than 20 sigma peak
brightprofpeak = np.max(baselineremoved[:,brightestprofile])
brightprofrms = rmsperepoch[brightestprofile]
if brightprofpeak/brightprofrms < 2 :
    resampled = scisig.resample(baselineremoved,bins/8)
    print 'resampling 8'
else:
    resampled = baselineremoved
baselineremoved = resampled

# Align data by the following way: x-correlate profiles with brightest, shift then x-correlate with average
aligned_data, template = Vf.aligndata(baselineremoved, brightestprofile, pulsar)
originaltemplate = np.copy(template)
suffix = ''
plt.plot(originaltemplate)
plt.savefig('./{0}/{0}_originaltemplate.png' .format(pulsar))
plt.clf()

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

# File that manually sets the start and end bin window of the pulse

f = open('/data1/pulsar_windows.txt','r')
pulsar_windows = []
pulsar_windows_lines = f.readlines()

for i in range(len(pulsar_windows_lines)):
    pulsar_windows.append(pulsar_windows_lines[i].split())
    if  pulsar_windows[-1][0] == pulsar:
        if len(pulsar_windows[-1]) == 3:
            binstartzoom = [int(pulsar_windows[-1][1])]
            binendzoom = [int(pulsar_windows[-1][2])]
            regioncounter = 1
        if len(pulsar_windows[-1]) == 5:
            binstartzoom = [int(pulsar_windows[-1][1]),int(pulsar_windows[-1][3])]
            binendzoom = [int(pulsar_windows[-1][2]),int(pulsar_windows[-1][4])]
            regioncounter = 2

print 'Found ', regioncounter, ' region(s)'
left = np.array(binstartzoom)
right = np.array(binendzoom)
binline = np.zeros(3)

# make sure they are in ascending order
binstartzoom.sort()
binendzoom.sort()

all_on_pulse = []
for b in range(regioncounter):
    all_on_pulse.append(range(binstartzoom[b],binendzoom[b]))

[item for sublist in all_on_pulse for item in sublist]

all_on_pulse = all_on_pulse[0]

on_pulse_data = aligned_data[all_on_pulse,]

# Normalise data to peak if flag -n used. Peak should already be aligned and at nbins/4:
if (args.normalise):
    aligned_data = Vf.norm_to_peak(aligned_data,on_pulse_data)
    suffix = '_norm'
    flagged = 0

    if not (os.path.exists('./{0}/flagged_profiles'.format(pulsar))):
        os.mkdir('./{0}/flagged_profiles'.format(pulsar))

    for i in range(aligned_data.shape[0]):
        originaltemplate[i]=np.median(aligned_data[i,:])
        
    med_peak = np.max(originaltemplate)
    originaltemplate = originaltemplate/med_peak
    aligned_data = aligned_data/med_peak
    # If profiles deviate far from the median profile, they are flagged and saved in a folder called 'flagged_profiles'
    sum_diff = np.zeros((aligned_data.shape[1]))

    for i in range(aligned_data.shape[1]):
        for k in range(aligned_data.shape[0]):
                sum_diff[i]+=abs(aligned_data[k,i]-originaltemplate[k])
            
    for j in range(aligned_data.shape[1]):
        if sum_diff[j] > 1.5*np.mean(sum_diff[:]):
            flagged+=1
            plt.plot(aligned_data[:,j],'b-')
            plt.plot(originaltemplate[:],'r-')
            plt.text(bins/3, 0.9,'Deviation from Median Profile: {0:.4f}' .format(sum_diff[j],fontsize=8))
            plt.text(bins/3, 0.8,'Mean of Profile Deviation: {0:.4f}' .format(np.mean(sum_diff[:]),fontsize=8))
            plt.text(bins/3, 0.7,'STDev of Profile Deviation: {0:.4f}' .format(np.std(sum_diff[:]),fontsize=8))
            plt.savefig('./{0}/flagged_profiles/{1}.png' .format(pulsar,mjdout[j],j))
            plt.clf()

    print "Percentage of profiles flagged as unusual:",np.round(float(flagged)/aligned_data.shape[1]*100.0,2),"% (",flagged,'out of',aligned_data.shape[1],")"

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

if (args.goodprofiles):
    dir='good_profiles{0}' .format(suffix)

    if regioncounter > 1:
        # Make interpulse plots of good profiles if needed
        if (args.normalise):
            Vf.goodplots_ip(pulsar,aligned_data[left[0]:right[0],:], aligned_data[left[1]:right[1],:],mjdout,dir,bins,left[1],originaltemplate[left[0]:right[0]], originaltemplate[left[1]:right[1]],-0.5*np.std(aligned_data[left[0]:right[0],:]),np.max(aligned_data[left[0]:right[0],:]),peakindex=bins/4-left[0])
        else: 
            Vf.goodplots_ip(pulsar,aligned_data[left[0]:right[0],:], aligned_data[left[1]:right[1],:],mjdout,dir,bins,left[1],originaltemplate[left[0]:right[0]], originaltemplate[left[1]:right[1]],-0.5*np.std(aligned_data[left[0]:right[0],:]),np.max(aligned_data[left[0]:right[0],:]),peakindex=bins/4-left[0],cal=1)
    else:
        # Make plots of good profiles if needed. cal = y or nothing
        if (args.normalise):
            Vf.makeplots(pulsar,aligned_data[left[0]:right[0],:],mjdout,dir,bins,template=originaltemplate[left[0]:right[0]],yllim=-0.1*np.max(originaltemplate),yulim=np.max(aligned_data[left[0]:right[0],:]),peakindex=bins/4-left[0])
        else:
            Vf.makeplots(pulsar,aligned_data[left[0]:right[0],:],mjdout,dir,bins,template=originaltemplate[left[0]:right[0]],yllim=-0.1*np.max(originaltemplate),yulim=np.max(aligned_data[left[0]:right[0],:]),peakindex=bins/4-left[0],cal=1)

# make file of median and std dev in each bin

if (args.normalise): 

    open('{0}/{0}_med.txt' .format(pulsar), 'w').close()
    open('{0}/{0}_std.txt' .format(pulsar), 'w').close()

    for i in range(aligned_data.shape[0]):
        bin_med = np.median(aligned_data[i,:])
        bin_std = np.std(aligned_data[i,:])
        f = open('{0}/{0}_med.txt' .format(pulsar), 'a')
        f.write('{0}\n' .format(bin_med))
        g = open('{0}/{0}_std.txt' .format(pulsar), 'a')
        g.write('{0}\n' .format(bin_std))
    f.close()
    g.close()

# Find the off-pulse standard deviation for use in Vgp.py

off_pulse_data = np.delete(aligned_data,all_on_pulse,0)

std_off_pulse_data = np.std(off_pulse_data)

z = open('{0}/{0}_off_pulse_std{1}.txt' .format(pulsar,suffix), 'w')
z.write('{0}' .format(std_off_pulse_data))

#    stdev = np.loadtxt('{0}/{0}_std.txt' .format(pulsar))

#    plt.plot(stdev)
#    plt.savefig('{0}_stdev.png' .format(pulsar))

# Make plots of removed profiles if needed
if (args.badprofiles):
    dir='removed_profiles{0}' .format(suffix)
    Vf.makeplots(pulsar,removedprofiles,mjdremoved,dir,bins,template=originaltemplate,yllim=-0.1*np.max(originaltemplate),yulim=np.max(aligned_data[left[0]:right[0],:]))

outputfile = '{0}/mjd{1}.txt' .format(pulsar,suffix)
np.savetxt(outputfile, mjdout)
outputfile = '{0}/mjdremoved{1}.txt' .format(pulsar,suffix)
np.savetxt(outputfile, mjdremoved)

# Save data from certain profiles for later plots

#for i in profs: 
#    print "I HERE",i
#    np.savetxt('./{0}/{0}_prof{1}.dat' .format(pulsar,i), aligned_data[:,i])
np.savetxt('./{0}/{0}_prof0.dat' .format(pulsar), aligned_data[:,0])


if (args.diagnosticplots):
    plt.imshow(data[:,:-1],cmap=plt.cm.binary,aspect = 'auto')
    plt.ylabel('Pulse Phase Bin')
    plt.xlabel('Observation Index')
#    plt.colorbar(orientation="horizontal")
    plt.savefig('./{0}/{0}_rawdata.png' .format(pulsar))
    plt.imshow(baselineremoved,cmap=plt.cm.binary,aspect = 'auto')
    plt.ylabel('Pulse Phase Bin')
    plt.xlabel('Observation Index')
    plt.savefig('./{0}/{0}_baselined.png' .format(pulsar))
    plt.clf()
